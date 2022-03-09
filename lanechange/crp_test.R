# This is an example where we run the Chinese Restaurant Process to perform clustering
source("~/Documents/anomaly-detection-ngsim/lanechange/crp_lib.R")
# ================
#
# GENERATE A SAMPLE DATASET
#
# =================
pi <- c(0.2, 0.6, 0.2)
mean <- c(0.3, 3.0, 9.0)
sd <- c(0.1, 0.2, 0.1)
n <- 100

x <- rep(0, n)
latent_cluster <- rep(0,n)

for (i in 1:n)
{
  cluster <- sample(1:3, replace=T, prob=pi)
  latent_cluster[i] <- cluster
  val <- rnorm(1, mean=mean[cluster], sd=sd[cluster])
  x[i] <- val
}

hist(x, breaks=50)

# ================
#
# DIRICHLET PROCESS GAUSSIAN MIXTURE MODEL (DP-GMM)
#
# =================

# Cluster base distribution
mu0 <- 1.0
var0 <- 4.0
tau0 <- solve(var0)

# Measurement variance
vary <- 0.05
tauy <- solve(vary)

# Attribute all datapoints with the first cluster
attr <- rep(1, n)
# Number of cluster
Nattr <- 1

# Cluster spawning coefficient
alpha <- 3

max_iter <- 1000 # max iter per clustering attempt

pb = txtProgressBar(min = 0, max = max_iter, initial = 0) 
for (iter in 1:max_iter)
{
  setTxtProgressBar(pb, iter)
  for (point_index in 1:length(x))
  {
    seatings <- n_customers_at_each_table(attr)
    x_of_customer <- x[point_index]
    table_of_current_customer <- attr[point_index]
    # First we un-seat the current customer (or segment) from his table
    seatings[table_of_current_customer] <- seatings[table_of_current_customer] - 1
    attr[point_index] <- -1
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
        nk <- seatings[table]
        tauk <- solve(current_var)
        # Sum of data in the current cluster
        cluster_sum <- sum(x[which(attr==table)])
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
        p_sit_at_table_k <- p_sit_at_occupied * dnorm(x_of_customer, mean=new_mean, sd=new_sd)
        p_sit_at_table[table] <- p_sit_at_table_k
      } else 
      {
        # Open a new table
        p_sit_at_unoccupied <- alpha / (N_seatings - 1 + alpha)
        new_sd <- sqrt(vary + var0)
        p_sit_at_table_Kplus1 <- p_sit_at_unoccupied * dnorm(x_of_customer, mean=mu0, sd=new_sd)
        p_sit_at_table[table] <- p_sit_at_table_Kplus1
      }
    } # end of loop through all tables
    # Normalize the p_sit_at_table
    p_sit_at_table <- p_sit_at_table / sum(p_sit_at_table)
    # Sample the new table for this customer to sit (including opening new one)
    table_to_sit <- sample(1:(Nattr+1), 1, replace=T, prob=p_sit_at_table)
    attr[point_index] <- table_to_sit
    if (table_to_sit == Nattr+1)
    {
      # The customer has decided upon opening a new table
      Nattr <- Nattr + 1 # increase the number of tables (i.e., clusters) by 1
    }
  }
  # cat("\r Iter ", iter)
}
close(pb)