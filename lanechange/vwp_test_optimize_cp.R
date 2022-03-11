# This is an example where we try to detect the "slopes"
# ================
#
# GENERATE A SAMPLE SIGNAL
#
# =================
source("~/Documents/anomaly-detection-ngsim/lanechange/vwp_lib.R")
x_sample <- rep(0,30)
x_mat <- matrix(, nrow=0, ncol = length(x_sample))
slopes <- c(1,3,1,3,1,3,1,3)
zero_period <- 12
rise_period <- 5
for (s in 1:20)
{
  slope_id <- 1 # to use different slopes at different rising periods
  for (t in 1:length(x_sample))
  {
    noise = rnorm(1,mean=0,sd=0.2)
    if (t %% zero_period > 0)
    {
      if (t %% 12 < rise_period + 1)
      {
        # this is a rising period
        x_sample[t] <- slopes[slope_id] * (t %% zero_period) %% rise_period + noise
      }
      else if (t %% 12 == rise_period + 1)
      {
        slope_id <- slope_id + 1
      }
      else {
        # this is a zero period
        x_sample[t] <- noise
      }
    }
  }
  x_mat <- rbind(x_mat, x_sample)
}

plot(x_sample)

# ================
#
# DEFINE A BASIS FUNCTION
#
# =================
basis = 1:30 
plot(basis)

# ================
#
# CALCULATE THE A (REPRESENTATION) AND HOMOGENEITY PROBA FOR
# ALL SEGMENTS IN EACH TIME SERIES (for quicker access)
#
# ================

# Measurement variance
vary <- 0.05
tauy <- solve(vary)
# To store all a values corresponding to segments beginning and terminating anywhere in each time series, for all time series
a_of_all_segments <- array(, c(ncol(x_mat), ncol(x_mat), nrow(x_mat)))
homo_of_all_segments <- array(, c(ncol(x_mat), ncol(x_mat), nrow(x_mat)))

for (ts_index in 1:nrow(x_mat))
{
  x_series <- x_mat[ts_index,]
  # We calculate the probability for each time series
  for (i in 1:ncol(x_mat))
  {
    for (j in i:ncol(x_mat))
    {
      # xij means xi and xj are included
      segment <- x_series[i:j]
      segment_a <- least_square(segment, basis)$b
      segment_homo <- log_p_homogeneity(segment, basis, segment_a, vary)
      a_of_all_segments[i,j,ts_index] <- segment_a
      homo_of_all_segments[i,j,ts_index] <- segment_homo
    }
  }
}


# ================
#
# DIRICHLET PROCESS GAUSSIAN MIXTURE MODEL (DP-GMM)
#
# =================

# Cluster base distribution
mu0 <- 1.0
var0 <- 1.0
tau0 <- solve(var0)

# All the changepoints position
cp <- list()
# Number of changepoints in each time series
Ncp <- get_number_of_changepoints(cp)
Ncp_max <- max(Ncp)
alpha_cp <- 2 # rate at which the NCP customer will open a new table

# For each time series, we will run a Chinese Restaurant Process to determine if we would like
# to add another changepoint
# The justificiation is for multiple time series, the number of changepoints should not vary so much
# between time series (is this reasoning good?)

for (ts in 1:nrow(x_mat))
{
  seatings_cp <- n_customers_at_each_table_cp(Ncp, 10)
  table_of_current_cp <- Ncp[ts]
  # We first unseat the NCP customer
  seatings_cp[table_of_current_cp] <- seatings_cp[table_of_current_cp] - 1
  N_seatings_cp <- sum(seatings_cp)
  for (table_cp in 1:(Ncp_max+1)) # table_cp: number of changepoints there are in 
  {
    nk <- seatings_cp[table_cp] # number of time series having the number of changepoint equals to table_cp
    if (table_cp <= Ncp_max)
    {
      # The NCP customer decides to sit with others
      p_ncp_sit_with_others <- nk / (N_seatings_cp - 1 + alpha_cp)
      
      # We calculate the search space of the changepoints placement
      scaling_coeff_for_nrow <- 1
      for (i in 1:table_cp)
      {
        scaling_coeff_for_nrow <- scaling_coeff_for_nrow * 1/i
      }
      total_num_of_nrow <- 1
      n_ts_length <- ncol(x_mat) # length of each time series
      for (i in 1:table_cp)
      {
        total_num_of_nrow <- total_num_of_nrow * (n_ts_length - i - 1)
      }
      # Begin filling in the proba for each changepoint placement hypothesis
      remaining_cps <- table_cp
      cp_placement <- list()
      
      while (remaining_cps > 0)
      {
        current_cp_index <- 1 + table_cp - remaining_cps
        current_cp_range <- cp_range(n_ts_length, current_cp_index, 0)
        for (i in current_cp_range)
        {
          
        }
      }
      
      
    }
  }
}

# We create an x_cut object which contains all the segments when cutting the x_mat time series
x_cut <- list()

# We create an "attribution" object which indicates the cluster that each segment belongs to
attr <- matrix(, nrow=nrow(x_mat), ncol=10) # 10 is the maximum number of changepoints in each time series

# Originally, this object will just be a duplicate of x_mat, all the original time series not cut
for (ts in 1:nrow(x_mat))
{
  x_cut <- append(x_cut, list(x_mat[ts,]))
}

# Then we place random changepoints among the time series
for (ts_index in 1:length(x_cut))
{
  ts_length <- length(x_cut[[ts_index]])
  change_at <- sample(2:ts_length-1, 1)
  # cp[ts_index] <- list(append(cp[[ts_index]], change_at))
  cp[ts_index] <- list(change_at)
}

# We can cut our original time series according to cp
xa_cut <- cut_x_with_cp(x_cut, cp, basis)
x_cut <- xa_cut$x_cut
a_cut <- xa_cut$a_cut

max_iter <- 1000 # max iter per clustering attempt
pb = txtProgressBar(min = 0, max = max_iter, initial = 0)
print('Running the Restaurant Process...')
for (iter in 1:max_iter)
{
  setTxtProgressBar(pb, iter)
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
          new_mean <- mu0
          new_sd <- sqrt(vary + var0)
          p_sit_at_table_Kplus1 <- p_sit_at_unoccupied * dnorm(a_of_customer, mean=mu0, sd=new_sd)
          p_sit_at_table[table] <- p_sit_at_table_Kplus1
        }
      } # end of loop through all tables
      # Normalize the p_sit_at_table
      p_sit_at_table <- p_sit_at_table / sum(p_sit_at_table)
      # Sample the new table for this customer to sit (including opening new one)
      table_to_sit <- sample(1:(Nattr+1), 1, replace=T, prob=p_sit_at_table)
      # Noting the new_mean of the table/cluster where our new customer will be sitting at
      cmu[table_to_sit] <- new_mean
      # Consider noting the variance also but since we have fixed variance, this should be fine
      attr[ts_index,segment_index] <- table_to_sit
      if (table_to_sit == Nattr+1)
      {
        # The customer has decided upon opening a new table
        Nattr <- Nattr + 1 # increase the number of tables (i.e., clusters) by 1
      }
    }
  } # end of CRP clustering
  
} # end of iterations
close(pb)