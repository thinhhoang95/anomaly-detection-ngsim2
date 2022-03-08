# This is an example where we try to detect the "slopes"
# ================
#
# GENERATE A SAMPLE SIGNAL
#
# =================
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

# Cluster base distribution
mu0 <- 3.0
var0 <- 3.0
tau0 <- solve(var0)

# Cluster params
cmu <- c(mu0)
cvar <- c(var0)

# All the changepoints
cp <- list()

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
# Assign the initial cluster of each segment, just assign to the same cluster for now
attr[, 1] <- 1
attr[, 2] <- 1
# Number of cluster
Nattr <- 1
# Cluster spawning coefficient
alpha <- 7

max_iter <- 15000 # max iter per clustering attempt

for (ts_index in 1:length(x_cut))
{
  for (segment_index in a_cut[[ts_index]])
  {
    segment <- a_cut[[ts_index]][[segment_index]]
    segment_a <- least_square(segment, basis)$b
    seatings <- n_customers_at_each_table(attr)
    table_of_current_customer <- attr[ts_index, segment_index]
    # First we un-seat the current customer (or segment) from his table
    seatings[table_of_current_customer] <- seatings[table_of_current_customer] - 1
    attr[ts_index][segment_index] <- -1
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
    
    for (table in 1:(Nattr+1))
    {
      if (table <= Nattr)
      {
        # Sit with others
        # Calculate posterior values
        # The current mean and variance
        current_mean <- cmu[table]
        current_var <- cvar[table]
        nk <- seatings[table]
        tauk <- solve(current_var)
        # Run a loop to calculate the sum of data of the current cluster
        cluster_sum <- 0
        for (ts_index in 1:length(x_cut))
        {
          for (segment_index in 1:length(x_cut[[ts_index]]))
          {
            if (attr[ts_index][segment_index] == table)
            {
              cluster_sum <- cluster_sum + as.vector(a_cut[[ts_index]][[segment_index]])
            }
          }
        } # Calculation of cluster_sum completed
        
      } else 
      {
        # Open a new table
      }
      
      # Dirichlet
      ni <- seatings[table]
      n <- sum(seatings)
      p_sit_with_others <- ni/(n - 1 + alpha)
    }
  }
}