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
# DIRICHLET PROCESS GAUSSIAN MIXTURE MODEL (DP-GMM)
#
# =================



# == Initialization of the Dirichlet Process ==
# ==
# Cluster base distribution
mu0 <- 1.0
var0 <- 1.0
tau0 <- solve(var0)

# Measurement variance
vary <- 0.05
tauy <- solve(vary)

# Cluster params
# cmu <- c(mu0) # if we want we could store the cluster mean in an auxiliary variable, but let's just leave it out for the moment
# cvar <- c(var0) # we don't actually use this because in the simplest case, we just assume a fixed variance for each cluster

# All the changepoints
cp <- list()
Ncp <- 1

# We create an x_cut object which contains all the segments when cutting the x_mat time series
x_cut <- list()

# Originally, this object will just be a duplicate of x_mat, all the original time series not cut
for (ts in 1:nrow(x_mat))
{
  x_cut <- append(x_cut, list(x_mat[ts,]))
}

# The number of maximum segments one time series can have
Ncp_max <- 10
# We create an "attribution" object which indicates the cluster that each segment belongs to (this is for the second CRP
# where we group the segments into the cluster
attr <- matrix(, nrow=nrow(x_mat), ncol=Ncp_max) # 10 is the maximum number of segments in each time series
attr[, 1] <- 1
Nattr_cp <- max(attr, na.rm=T)
alpha_cp <- 2 # the rate at which the time series will decide to open a larger number of changepoints that other time series has not yet considered

# == End of Initialization of the Dirichlet Process ==
# ==


Ncp_of_each_ts <- get_number_of_changepoints(cp)

for (ts_index in 1:length(x_cut))
{
  # == Step 1: Setting a new number of changepoints for the current TS ==
  # Get the seatings of each number of changepoints
  cp_seatings <- n_customers_at_each_table_cp(Ncp_of_each_ts, Ncp_max)
  # We first unseat the current time series
  table_of_current_ts <- Ncp_of_each_ts[ts_index]
  cp_seatings[table_of_current_ts] <- cp_seatings[table_of_current_ts] - 1
  # We run a fake Chinese Restaurant Process for the current time series to choose
  # the number of changepoints cp_table
  for (cp_table in 1:(Nattr_cp+1))
  {
    if (cp_table <= Nattr_cp)
    {
      # Calculate proba that the time series will open a higher number of changepoints that
      # other time series yet to consider
      nk_cp <- cp_seatings[cp_table]
      N_seatings_cp <- sum(cp_seatings)
      p_sit_with_others_cp <- nk_cp / (N_seatings_cp - 1 + alpha_cp)
    } else {
      # Calculate proba that the time series will open a higher number of changepoints that
      # other time series yet to consider
    }
  }
  # == End of step 1: Setting a new number of changepoints for the current TS ==
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
# Assign the initial cluster of each segment, just assign to the same cluster for now
attr[, 1] <- 1
attr[, 2] <- 1
# Number of cluster
Nattr <- 1
# Cluster spawning coefficient
alpha <- 1

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
}
close(pb)