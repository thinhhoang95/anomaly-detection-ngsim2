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

# We create an x_cut object which contains all the segments when cutting the x_mat time series
x_cut <- list()

for (ts in 1:nrow(x_mat))
{
  x_cut <- append(x_cut, list(x_mat[ts,]))
}

rows_to_skip_in_x_cut <- nrow(x_mat)
# We first consider a random, unoptimized changepoint placement at each time series
for (its in 1:length(x_cut))
{
  ts <- x_cut[[its]]
  ts_length = length(ts)
  change_at <- sample(1:ts_length, 1)
  if (change_at < 2 | change_at > ts_length - 2)
  {
    # the sampled changepoint is too close to the beginning or the end of the time series
    # just leave the time series as is
    x_cut <- append(x_cut, list(ts))
  } 
  else
  {
    x_cut <- append(x_cut, list(ts[1:(change_at - 1)]))
    x_cut <- append(x_cut, list(ts[change_at:ts_length]))
  }
}

# Since all the time series in the beginning has been cut, let's remove them from the 
# original x_cut
x_cut <- tail(x_cut, -rows_to_skip_in_x_cut)

# So up to this moment, we have a bunch of segments randomly cut from the dataset,
# we will calculate its representations (a*)

# 
a_cut <- list()
a_likely_cut <- list()
assignment <- rep(1,length(x_cut))
cluster_means <- list(c(0.))
cluster_var <- list(c(0.2))

# We can start the Chinese Restaurant Process to begin the clustering process


tau_x <- 10 # 1/tau_x = sigma_x^2, precision of measurement - or model's belief about the noise

for (i in 1:length(x_cut))
{
  ts <- x_cut[[i]]
  ts_length = length(ts)
  lss = least_square(ts, basis)
  # Regression coefficients
  lssb <- as.vector(lss$b) # a*
  a_cut <- append(a_cut, list(lssb))
  # Homogeneity measurement: likelihood of a* generating measurements
  # First, reconstruct the "supposed" ts from a*
  ts_ref <- as.vector(reconstruct(ts, lssb, basis))
  # Compare the ts and ts_ref time by time
  a_likely_cut[i] <- 1./ts_length * (-0.5) * tau_x * sum((ts_ref - ts)^2)
}
