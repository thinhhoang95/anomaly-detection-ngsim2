# =======================
#
# DATA PREPARATION
#
# =======================
# Load CSV that contains lane changing data
data <- as.matrix(unname(read.csv("lanechange.csv", header=F)))

# We try to resample all the time series to reduce the search space
data_rs <- data[,seq(1,ncol(data), 20)]

# Trajectory 101 includes a lane change
plot(data[101,], type = 'l')
plot(data_rs[101,], type = 'l', col='blue')

# Select a subset of trajectories for segmentation learning
selected_rows <- seq(1,nrow(data),3)
cat("Selected rows: ", selected_rows, "\n")
data_srs <- data_rs[selected_rows,]

cat("Dataset dimensions: ", dim(data_srs))

# Plot all the trajectories to visualize first
for (i in 1:nrow(data_srs))
{
  # cat("Plotting graph ", i, "\n")
  if (i==1)
  {
    plot(data_srs[1,], ylim=range(data_srs), type='l')
  } else {
    lines(data_srs[i,])
  }
}

x_mat <- data_srs # rename the dataset so that it is compatible with the code we wrote earlier

# ================
#
# DEFINE A BASIS FUNCTION
#
# =================
basis = 0:30 

# ================
#
# LOAD THE VIETNAMESE WEDDING PROCESS LIBRARY
#
# =================
source("~/Documents/anomaly-detection-ngsim/lanechange/vwp_lib.R")

# ================
#
# DIRICHLET PROCESS GAUSSIAN MIXTURE MODEL (DP-GMM)
#
# =================



# == Initialization of the Dirichlet Process ==
# ==

# Cluster base distribution
mu0 <- 1.0
var0 <- 0.5
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
alpha <- 1e-5

# Cluster params
cmu <- rep(mu0, Ncluster_max)
cvar <- var0 # we don't actually use Inverse Wishart Prior because in the simplest case, we just assume a fixed variance for each cluster
ctau <- solve(cvar)

# We create an "attribution" object which indicates the cluster that each segment belongs to (this is for the second CRP
# where we group the segments into the cluster
attr <- matrix(, nrow=nrow(x_mat), ncol=Ncp_max) # 10 is the maximum number of segments in each time series
# Assign the initial cluster for each segment, just assign to the same cluster for now
attr[, 1] <- 1
attr[, 2] <- 1

# For the number of changepoints, Nattr_cp is the maximum number of changepoints (tables).
# Each time, a time series will be asked if they want to increase this number (open a new table).
Nattr_cp <- max(attr, na.rm=T)
alpha_cp <- 1e-1 # the rate at which the time series will decide to open a larger number of changepoints that other time series has not yet considered


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

# == End of Initialization of the Dirichlet Process ==
# ==

# == Begin the iterative procedure ==
# == >> << == 

max_big_iter <- 1
for (big_iter in 1:max_big_iter)
{
  all_pah <- list()
  cat('Iteration ', big_iter, '\n')
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
    # We first unseat the current time series
    table_of_current_ts <- Ncp_of_each_ts[ts_index]
    cp_seatings[table_of_current_ts] <- cp_seatings[table_of_current_ts] - 1
    # We run a fake Chinese Restaurant Process for the current time series to choose
    # the number of changepoints cp_table
    ts_length <- length(x_mat[ts_index,])
    pah <- precalculate_a_and_homogeneity(x_mat[ts_index,], basis, vary)
    all_pah[ts_index] <- list(pah) # We store all precalculated values so that the clustering can re use these values
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
        N_seatings_cp <- sum(cp_seatings)
        p_sit_with_others_cp <- nk_cp / (N_seatings_cp - 1 + alpha_cp)
        # We find the optimal segmentation given the number of changepoints equal cp_table
        part <- find_most_likely_partitioning_of_a_time_series(ts_length, cp_table, cmu, ctau, Nattr, pah)
        p_this_cp_is_correct <- exp(part$lp) # given the maximum number of changepoints and the cluster means
        # p_table_cp[cp_table] <- p_this_cp_is_correct * p_sit_with_others_cp
        # We use the log directly for better numerical stability
        p_table_cp[cp_table] <- part$lp + log(p_sit_with_others_cp)
        # cat("\n Likelihood of CRP: ", p_sit_with_others_cp)
        # cat("\n Likelihood of CP: ", part$lp)
        # cat("\n Proposed CP: ", part$cp)
        proposed_cp <- append(proposed_cp, list(part$cp))
        proposed_cluster_attr <- append(proposed_cluster_attr, list(part$attr))
      } else {
        # Calculate proba that the time series will open a higher number of changepoints that
        # other time series yet to consider
        p_open_a_new_table_cp <- alpha_cp / (N_seatings_cp - 1 + alpha_cp)
        part <- find_most_likely_partitioning_of_a_time_series(ts_length, cp_table, cmu, ctau, Nattr, pah)
        p_this_cp_is_correct <- exp(part$lp) # given the maximum number of changepoints and the cluster means
        # p_table_cp[cp_table] <- p_this_cp_is_correct * p_open_a_new_table_cp
        # We use the log directly for better numerical stability
        p_table_cp[cp_table] <- part$lp + log(p_open_a_new_table_cp)
        # cat("\n Likelihood of CRP: ", p_open_a_new_table_cp)
        # cat("\n Likelihood of CP: ", part$lp)
        # cat("\n Proposed CP: ", part$cp)
        proposed_cp <- append(proposed_cp, list(part$cp))
        proposed_cluster_attr <- append(proposed_cluster_attr, list(part$attr))
      }
    } # end of looping through all the tables
    # We scale the largest value to 1
    p_table_cp <- p_table_cp + max(p_table_cp)
    # print(p_table_cp)
    p_table_cp <- exp(p_table_cp)
    # Then we normalize this vector, which is always possible since there guarantees
    # to exist a number 1 in the p_table_cp vector, so it can't divide by 0
    p_table_cp <- p_table_cp / sum(p_table_cp)
    
    #cat("\n Table assignment proba: ", p_table_cp)
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
    best_cluster_attr <- proposed_cluster_attr[[table_to_sit_cp]]
    attr[ts_index, 1:length(best_cluster_attr)] <- best_cluster_attr
    
    # == End of step 1: Setting a new number of changepoints for the current TS ==
  }
  close(pb)
  
  cat("Maximum changepoints is now ", Nattr_cp, "\n")
  
  # == Step 2: Assigning cluster for each segment ==
  # NOTICE: Consider integrating this loop into Step 1's loop to have better handling of memory
  
  # We can cut our original time series according to cp
  xa_cut <- cut_x(x_mat, cp, basis, all_pah)
  x_cut <- xa_cut$x_cut
  a_cut <- xa_cut$a_cut
  
  max_iter <- 700 # max iter per clustering attempt
  pb2 = txtProgressBar(min = 0, max = max_iter, initial = 0)
  print('Clustering segments...')
  for (iter in 1:max_iter)
  {
    setTxtProgressBar(pb2, iter)
    for (ts_index in 1:length(x_cut))
    {
      for (segment_index in 1:length(a_cut[[ts_index]]))
      {
        seatings <- n_customers_at_each_table(attr)
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
            for (aux_ts_index in 1:length(x_cut))
            {
              for (aux_segment_index in 1:length(x_cut[[aux_ts_index]]))
              {
                if (attr[aux_ts_index,aux_segment_index] == table)
                {
                  cluster_sum <- cluster_sum + as.vector(a_cut[[aux_ts_index]][[aux_segment_index]])
                }
              }
            } # Calculation of cluster_sum completed
            taup <- nk * tauk + tau0
            sigp <- solve(taup)
            new_mean <- sigp * (tauk * cluster_sum + mu0 * tau0)
            cmu[table] <- new_mean
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
} # end of big iter