# VIETNAMESE WEDDING PROCESS FOR TRAJECTORIES OF A PARTICULAR EDGE
# Data is generated from trajectory_preprocessor.ipynb

# Housekeeping
rm(list = ls(all.names = TRUE))
gc()

df <- read.csv('edge_trajectories.csv')
cat('CSV file loaded \n')
# For each vehicle in the dataset, we will generate a matrix
vehicle_ids <- unique(df$vehicle_id)
cat(length(vehicle_ids), " vehicles were found in the dataset \n")
result <- matrix(NA, nrow=length(vehicle_ids), ncol=99) # 99 comes from the Python notebook
i <- 1
for (vehicle_id in vehicle_ids)
{
  row <- df[df$vehicle_id==vehicle_id,]$x
  result[i,] <- tail(row, 1)
  result[i,1:length(row)] <- row
  i <- i+1
}

result <- result[1:30,] # just take the first 30 rows

# Plot the trajectories
for (i in 1:nrow(result))
{
  if (i==1)
  {
    plot(1:99, result[i,], type='l', xlab = 'Timestep', ylab = 'Coordinate (m)')
  } else {
    lines(1:99, result[i,])
  }
}
# Now the matrix result contains all trajectories of the specified edge
# Next we run the VWP process to segment out the dominant modes in the trajectories, then later,
# FPCA is used to obtain the best linear basis to learn the representation's distribution

cat('Dataset dimensions: ', dim(result), '\n')

# This example uses multivariate Student t distribution rather than fixed variance of each cluster
# this would be a more rigorous approach than the fixed variance one (although the fixed variance is much simpler)
# and equivalently effective

source("vwp_visualize.R")
debugSource("vwp_libc_wishart_wishart.R")
library(mvtnorm) # multivariate Student t distribution: for dmvt
# source("~/Documents/anomaly-detection-ngsim/lanechange/vwp_lib.R")

x_mat <- result
rm(result) #house-keeping

# ================
#
# DEFINE A BASIS FUNCTION
#
# =================
basis = 0:101

# ================
#
# DIRICHLET PROCESS GAUSSIAN MIXTURE MODEL (DP-GMM)
#
# =================



# == Initialization of the Dirichlet Process ==
# ==

# Cluster base distribution
mu0 <- 0.0
var0 <- 0.1 # cluster "concentration"
tau0 <- solve(var0)

# Measurement variance
vary <- 0.05
tauy <- solve(vary)

# All the changepoints
cp <- list()

# We create an x_cut object which contains all the segments when cutting the x_mat time series
x_cut <- list()

# Originally, this object will just be a duplicate of x_mat, all the original time series not cut
for (ts in 1:nrow(x_mat))
{
  x_cut <- append(x_cut, list(x_mat[ts,]))
}

# The number of maximum segments one time series can have
Ncp_max <- 10
# The maximum number of clusters each segment can take
Ncluster_max <- 10

# Number of current clusters (i.e., tables that all segments might sit)
Nattr <- 1
# Cluster spawning coefficient
alpha <- 1e-18

# Cluster params
cmu <- rep(mu0, Ncluster_max)
cvar <- rep(var0, Ncluster_max) # we don't actually use Inverse Wishart Prior because in the simplest case, we just assume a fixed variance for each cluster
ctau <- 1/cvar # the precision for each cluster

# We create an "attribution" object which indicates the cluster that each segment belongs to (this is for the second CRP
# where we group the segments into the cluster
attr <- matrix(1, nrow=nrow(x_mat), ncol=Ncp_max) # 10 is the maximum number of segments in each time series
# (We also assign cluster 1 to all the segments in the beginning)
# Assign the initial cluster for each segment, just assign to the same cluster for now
# attr[, 1] <- 1
# attr[, 2] <- 1

# For the number of changepoints, Nattr_cp is the maximum number of changepoints (tables).
# Each time, a time series will be asked if they want to increase this number (open a new table).
Nattr_cp <- max(attr, na.rm=T)
alpha_cp <- 1e-256 # the rate at which the time series will decide to open a larger number of changepoints that other time series has not yet considered


# We place random changepoints among the time series
# the actual positions of the changepoints are not important
# because we will change them anyway. What matters is many time series will originally sit at
# the same table. This influences the probability to open new table.
for (ts_index in 1:length(x_cut))
{
  ts_length <- length(x_cut[[ts_index]])
  change_at <- sample(2:ts_length-1, 1)
  # cp[ts_index] <- list(append(cp[[ts_index]], change_at))
  cp[ts_index] <- list(change_at)
}

# For dynamic programming of changepoint positions, this is the blank solution to start
# the algorithm will gradually write into the three arrays contained in this variable (default_sol)
default_sol <- list("p" = array(NA, dim=c(ts_length, Ncp_max)), "cp" = array(NA, dim=c(Ncp_max, ts_length, Ncp_max)), "a" = array(NA, dim=c(Ncp_max, ts_length, Ncp_max)))

# == End of Initialization of the Dirichlet Process ==
# ==

# == Begin the iterative procedure ==
# == >> << == 

# Precalculate the homogeneity and regression slopes
pa_homo <- precalculate_pa_homo(x_mat, basis, vary)

# We obtain the first x_cut and a_cut
# We can cut our original time series according to cp
xa_cut <- cut_x(x_mat, cp, basis, pa_homo)
x_cut <- xa_cut$x_cut
a_cut <- xa_cut$a_cut
l_cut <- xa_cut$l_cut

# Normal Inverse Wishart Hyperparameters
# For the slopes
hypa <- list("mu" = mu0, "nu" = 1, "kappa" = 5, "lambda" = 0.5)
# For the segment length
hypa_l <- list("mu" = 20, "nu" = 1, "kappa" = 2, "lambda" = 0.5)

max_big_iter <- 1
for (big_iter in 1:max_big_iter)
{
  all_pah <- list()
  cat('\n Iteration ', big_iter, '\n')
  print('Optimizing changepoints...')
  pb = txtProgressBar(min = 0, max = dim(x_mat)[1], initial = 0)
  for (ts_index in 1:(dim(x_mat)[1]))
  {
    setTxtProgressBar(pb, ts_index)
    # cat("\n ==========")
    # cat("\n Processing time series ", ts_index)
    
    # == Step 1: Setting a new number of changepoints for the current TS ==
    
    # Get the seatings of each number of changepoints
    Ncp_of_each_ts <- get_number_of_changepoints(cp)
    cp_seatings <- n_customers_at_each_table_cp(Ncp_of_each_ts, Ncp_max)
    N_seatings_cp <- sum(cp_seatings)
    # We first unseat the current time series
    table_of_current_ts <- Ncp_of_each_ts[ts_index]
    cp_seatings[table_of_current_ts] <- cp_seatings[table_of_current_ts] - 1
    # We run a fake Chinese Restaurant Process for the current time series to choose
    # the number of changepoints cp_table
    ts_length <- length(x_mat[ts_index,])
    # pah <- precalculate_a_and_homogeneity(x_mat[ts_index,], basis, vary)
    # all_pah[ts_index] <- list(pah) # We store all precalculated values so that the clustering can re use these values
    p_table_cp <- rep(0, Nattr_cp+1)
    proposed_cp <- list()
    proposed_cluster_attr <- list()
    
    for (cp_table in 1:(Nattr_cp+1))
    {
      if (cp_table <= Nattr_cp)
      {
        # Calculate proba that the time series will open a higher number of changepoints that
        # other time series yet to consider
        nk_cp <- cp_seatings[cp_table]
        p_sit_with_others_cp <- nk_cp / (N_seatings_cp - 1 + alpha_cp)
        if (p_sit_with_others_cp <= 0) # this corresponds to no changepoint was detected
        {
          p_sit_with_others_cp = 0 # taking the log of this will yield -Inf
        } 
        
        # We find the optimal segmentation given the number of changepoints equal cp_table
        segment_seatings <- n_customers_at_each_table(attr, Nattr)
        # Get the clusters' t-distribution params
        clusters <- get_clusters_to_use_in_dynamic_programming(Nattr, attr, a_cut, l_cut, hypa, hypa_l, ts_index)
        dp_sol <- dynamic_program_cp(default_sol, cp_table, x_mat[ts_index,], ts_index, pa_homo, 1, basis, vary, Nattr, clusters, segment_seatings, alpha)
        dp_cp <- dp_sol$cp[1:cp_table,1,cp_table]
        dp_lp <- dp_sol$p[1,cp_table]
        
        p_this_cp_is_correct <- exp(dp_lp) # given the maximum number of changepoints and the cluster means
        # p_table_cp[cp_table] <- p_this_cp_is_correct * p_sit_with_others_cp
        # We use the log directly for better numerical stability
        p_table_cp[cp_table] <- dp_lp + log(p_sit_with_others_cp)
        
        if (any(is.na(p_table_cp)))
        {
          print("NA detected")
        }
        
        # cat("\n NCP = ", cp_table, "; CP: ", dp_cp)
        # cat("\n Likelihood of CRP: ", p_sit_with_others_cp)
        # cat("\n Likelihood of CP: ", dp_lp)
        
        proposed_cp <- append(proposed_cp, list(dp_cp))
        # proposed_cluster_attr <- append(proposed_cluster_attr, list(part$attr))
      } else {
        
        # Calculate proba that the time series will open a higher number of changepoints that
        # other time series yet to consider
        p_open_a_new_table_cp <- alpha_cp / (N_seatings_cp - 1 + alpha_cp)
        if (p_open_a_new_table_cp <= 0)
        {
          p_open_a_new_table_cp <- 1
        }
        # We find the optimal segmentation given the number of changepoints equal cp_table
        segment_seatings <- n_customers_at_each_table(attr, Nattr)
        # Get the clusters' t-distribution params
        clusters <- get_clusters_to_use_in_dynamic_programming(Nattr, attr, a_cut, l_cut, hypa, hypa_l, ts_index)
        dp_sol <- dynamic_program_cp(default_sol, cp_table, x_mat[ts_index,], ts_index, pa_homo, 1, basis, vary, Nattr, clusters, segment_seatings, alpha)
        dp_cp <- dp_sol$cp[1:cp_table,1,cp_table]
        dp_lp <- dp_sol$p[1,cp_table]
        p_this_cp_is_correct <- exp(dp_lp) # given the maximum number of changepoints and the cluster means
        # p_table_cp[cp_table] <- p_this_cp_is_correct * p_sit_with_others_cp
        # We use the log directly for better numerical stability
        p_table_cp[cp_table] <- dp_lp + log(p_open_a_new_table_cp)
        
        if (any(is.na(p_table_cp)))
        {
          print("NA detected")
        }
        
        # cat("\n *NCP = ", cp_table, "; CP: ", dp_cp)
        # cat("\n Likelihood of CRP: ", p_open_a_new_table_cp)
        # cat("\n Likelihood of CP: ", dp_lp)
        
        proposed_cp <- append(proposed_cp, list(dp_cp))
      }
      
    } # end of looping through all the tables
    # We scale the largest value to 1
    if (any(is.na(p_table_cp)))
    {
      print("NA detected")
    }
    p_table_cp <- p_table_cp - max(p_table_cp)
    # cat("\n LOG assignment proba: ", p_table_cp)
    # print(p_table_cp)
    p_table_cp <- exp(p_table_cp)
    # Then we normalize this vector, which is always possible since there guarantees
    # to exist a number 1 in the p_table_cp vector, so it can't divide by 0
    p_table_cp <- p_table_cp / sum(p_table_cp)
    if (any(is.na(p_table_cp)))
    {
      print("NA detected")
    }
    # cat("\n Table assignment proba: ", p_table_cp)
    # Sample the new table for this customer to sit (including opening new one)
    table_to_sit_cp <- sample(1:(Nattr_cp+1), 1, replace=T, prob=p_table_cp)
    cp_seatings[table_to_sit_cp] <- cp_seatings[table_to_sit_cp] + 1
    #cat("\n Choosing the number of changepoints as ", table_to_sit_cp)
    cp[ts_index] <- list(proposed_cp[[table_to_sit_cp]])
    if (table_to_sit_cp == (Nattr_cp + 1))
    {
      Nattr_cp <- Nattr_cp + 1
      # cat("\n Opened a new number of changepoints. Maximum changepoints is now ", Nattr_cp)
    } else {
      # cat("\n Keep the same number of changepoints. Maximum changepoints is now ", Nattr_cp)
    }
    # Assigning the cluster attribution
    # best_cluster_attr <- proposed_cluster_attr[[table_to_sit_cp]]
    # attr[ts_index, 1:length(best_cluster_attr)] <- best_cluster_attr
    
    # == End of step 1: Setting a new number of changepoints for the current TS ==
  }
  close(pb)
  
  cat("Maximum changepoints is now ", Nattr_cp, "\n")
  
  # == Step 2: Assigning cluster for each segment ==
  # NOTICE: Consider integrating this loop into Step 1's loop to have better handling of memory
  
  # We can cut our original time series according to cp
  xa_cut <- cut_x(x_mat, cp, basis, pa_homo)
  x_cut <- xa_cut$x_cut
  a_cut <- xa_cut$a_cut
  l_cut <- xa_cut$l_cut
  
  max_iter <- 500 # max iter per clustering attempt
  pb2 = txtProgressBar(min = 0, max = max_iter, initial = 0)
  print('Clustering segments...')
  # We start fresh, since we have a new batch of customers (as the changepoints changed, the old assignments became invalid)
  # we shall run the CRP to cluster all of these segments, and get the new means for the clusters
  cmu <- rep(mu0, Ncluster_max)
  cmu_l <- rep(0, Ncluster_max)
  Nattr <- 1
  attr <- matrix(1, nrow=nrow(x_mat), ncol=Ncp_max) # 10 is the maximum number of segments in each time series
  for (iter in 1:max_iter)
  {
    setTxtProgressBar(pb2, iter)
    for (ts_index in 1:length(x_cut))
    {
      for (segment_index in 1:length(a_cut[[ts_index]]))
      {
        seatings <- n_customers_at_each_table(attr, Nattr)
        a_of_customer <- a_cut[[ts_index]][[segment_index]]
        table_of_current_customer <- attr[ts_index, segment_index]
        # First we un-seat the current customer (or segment) from his table
        seatings[table_of_current_customer] <- seatings[table_of_current_customer] - 1
        attr[ts_index,segment_index] <- -1
        if (seatings[table_of_current_customer] == 0)
        {
          # If the unseated customer turns out to be the last one sitting at that table, we remove the table via a two steps process
          # First, we swap the (just emptied) table with the last table (i.e., table Nattr) by reassigning the seats of all customers sitting there
          # Then we let Nattr - 1 so the last empty table is gone
          attr[attr==Nattr] <- table_of_current_customer # resitting all the customers at the last table
          seatings[table_of_current_customer] <- seatings[Nattr] # swap the number of customers of the two tables being swapped
          Nattr <- Nattr - 1 # remove the last (now empty) table
        }
        # Now run the Chinese Restaurant Process
        # Each cluster will try to "engulf" the segment's a and see its mean and variance shifted
        # to posterior mean and variance values
        
        # This variable holds the proba of all tables that the customer is going to sit in
        p_sit_at_table <- rep(0, Nattr+1)
        for (table in 1:(Nattr+1))
        {
          if (table <= Nattr)
          {
            # Sit with others
            # Calculate posterior values
            # The current mean and variance
            # current_mean <- mu0 # actually not used
            current_var <- var0 # use the fixed variance for cluster in this version of code
            # a posterior with inverse-Wishart prior and likelihood from the data seems like a better choice
            nk <- seatings[table]
            tauk <- solve(current_var)
            # Run a loop to calculate the sum of data of the current cluster
            cluster_sum <- 0
            cluster_sum_l <- 0
            for (aux_ts_index in 1:length(x_cut))
            {
              for (aux_segment_index in 1:length(x_cut[[aux_ts_index]]))
              {
                if (attr[aux_ts_index,aux_segment_index] == table)
                {
                  cluster_sum <- cluster_sum + as.vector(a_cut[[aux_ts_index]][[aux_segment_index]])
                  cluster_sum_l <- cluster_sum_l + as.vector(l_cut[[aux_ts_index]][[aux_segment_index]])
                }
              }
            } # Calculation of cluster_sum completed
            taup <- nk * tauk + tau0
            sigp <- solve(taup)
            new_mean <- sigp * (tauk * cluster_sum + mu0 * tau0)
            cmu[table] <- new_mean
            cmu_l[table] <- cluster_sum_l / nk
            new_var <- sigp + vary
            new_sd <- sqrt(new_var)
            N_seatings <- sum(seatings)
            # The likelihood of "sitting at an occupied table" is
            p_sit_at_occupied <- nk / (N_seatings - 1 + alpha)
            # Replace a_of_customer is the line below if we want to also account for the homogeneity of the measurements
            # here it works like the most basic HDP-GMM (that doesn't even account for )
            p_sit_at_table_k <- p_sit_at_occupied * dnorm(a_of_customer, mean=new_mean, sd=new_sd)
            p_sit_at_table[table] <- p_sit_at_table_k
          } else 
          {
            # Open a new table
            p_sit_at_unoccupied <- alpha / (N_seatings - 1 + alpha)
            new_sd <- sqrt(vary + var0)
            p_sit_at_table_Kplus1 <- p_sit_at_unoccupied * dnorm(a_of_customer, mean=mu0, sd=new_sd)
            p_sit_at_table[table] <- p_sit_at_table_Kplus1
          }
        } # end of loop through all tables
        # Normalize the p_sit_at_table
        p_sit_at_table <- p_sit_at_table / sum(p_sit_at_table)
        # Sample the new table for this customer to sit (including opening new one)
        table_to_sit <- sample(1:(Nattr+1), 1, replace=T, prob=p_sit_at_table)
        attr[ts_index,segment_index] <- table_to_sit
        if (table_to_sit == Nattr+1)
        {
          # The customer has decided upon opening a new table
          Nattr <- Nattr + 1 # increase the number of tables (i.e., clusters) by 1
        }
      }
    }
  } # clustering loop
  close(pb2)
  cat("Maximum clusters is now ", Nattr, "\n")
  # We will try to merge changepoints that are not changepoints at all!
  # Because it assigns the segment after the changepoint to the same cluster
  # as the current segment
  print('Merging segments...')
  merge_result <- merge_segments(cp, attr)
  cp <- merge_result$cp
  attr <- merge_result$attr
  # Recalculate the maximum changepoints opened
  Nattr_cp <- get_Nattr_cp_from_attr(attr)
  cat('Maximum changepoints is now ', Nattr_cp, '\n')
  cat('Cluster means: ', cmu[1:Nattr], '\n')
  cat('Average length of cluster segments: ', cmu_l[1:Nattr], '\n')
} # end of big iter

# Graph drawing scripts

# Some trajectory data
for (i in 1:nrow(x_mat))
{
  if (i==1)
  {
    plot(1:ncol(x_mat), x_mat[i,], type='l', ylim=c(min(x_mat), max(x_mat)), xlab='Time sample (s)', ylab='Lateral coordinate (ft)')
  }
  lines(1:ncol(x_mat), x_mat[i,])
}

xa_cut <- cut_x(x_mat, cp, basis, pa_homo)
x_cut <- xa_cut$x_cut
a_cut <- xa_cut$a_cut
l_cut <- xa_cut$l_cut
# histogram of slopes
hist(unlist(a_cut), breaks=200, main="", xlab="Value", ylab="Frequency")

###
# ============================= # 
# Learn the decomposition basis #
# ============================= #
###

cluster_trajs <- matrix(NA, nrow = 0, ncol = 50)
i <- 1
for (ts_index in 1:length(x_cut))
{
  for (segment_index in 1:length(x_cut[[ts_index]]))
  {
    if (attr[ts_index,segment_index] == 2)
    {
      row <- rep(NA, 50)
      # Regress X-Y
      traj <- x_cut[[ts_index]][[segment_index]]
      traj <- traj - traj[1]
      t_traj <- 1:(length(traj))
      # Extrapolate the time series to have 50 samples
      model <- lm(traj ~ t_traj)
      dtraj <- data.frame(t_traj= (length(traj)+1):50)
      lengthened_traj <- predict(model, dtraj)
      row <- c(traj, as.vector(lengthened_traj))
      cluster_trajs <- rbind(cluster_trajs, row)
      i <- i + 1
    }
  }
}

# Plot the trajectories in the cluster
# Plot the trajectories
for (i in 1:nrow(cluster_trajs))
{
  if (i==1)
  {
    plot(1:50, cluster_trajs[i,], type='l', xlab = 'Timestep', ylab = 'Coordinate (m)')
  } else {
    lines(1:50, cluster_trajs[i,])
  }
}

# Perform FPCA decomposition
# ===
# 1. Create a tibble in the correct format for FPCA
t_vec <- vector(mode = "list", length = nrow(cluster_trajs))
x_vec <- vector(mode = "list", length = nrow(cluster_trajs))
for (i in 1:nrow(cluster_trajs))
{
  t_vec[i] <- list(1:50)
  x_vec[i] <- list(cluster_trajs[i,])
}
tibb <- tibble(t_vec, x_vec)
# 2. Perform FPCA decomposition
library(fdapace)
CheckData(tibb$x_vec,tibb$t_vec)
W_fpca <- FPCA(tibb$x_vec, tibb$t_vec)