visualize_belief <- function(x_mat, cp, cmu, attr, index)
{
  ts_length <- ncol(x_mat)
  cp_of_current_ts <- cp[[index]]
  cp_of_current_ts_extended <- c(1, cp_of_current_ts, ts_length)
  x_recons <- rep(0, ts_length)
  for (i in 1:(length(cp_of_current_ts)+1))
  {
    segment_starts_at <- cp_of_current_ts_extended[i]
    segment_ends_at <- cp_of_current_ts_extended[i+1]
    # cat("\n Segment starts at", segment_starts_at)
    # cat("\n Segment ends at", segment_ends_at)
    # Reconstruct the segment
    a <- cmu[attr[index, i]]
    for (t in segment_starts_at:segment_ends_at)
    {
      if (t==1)
      {
        x_recons[t] <- x_mat[index, t]
      } else {
        if (t==segment_starts_at + 1)
        {
          x_recons[t] <- x_mat[index, t]
        } else {
          x_recons[t] <- x_recons[t-1] + a
        }
      }
      
    }
  }
  
  plot(x_recons, col='red', pch=19, ylim=range(cbind(x_mat[index,], x_recons))+c(-2,2), xlab="Time (s)", ylab="Coordinate (ft)")
  points(x_mat[index,], col='green')
}

explain_belief <- function(x_mat, cp, cmu, attr, index)
{
  visualize_belief(x_mat, cp, cmu, attr, index)
  ts_length <- ncol(x_mat)
  cp_of_current_ts <- cp[[index]]
  cp_of_current_ts_extended <- c(1, cp_of_current_ts, ts_length + 1)
  x_recons <- rep(0, ts_length)
  for (i in 1:(length(cp_of_current_ts)+1))
  {
    segment_starts_at <- cp_of_current_ts_extended[i]
    segment_ends_at <- cp_of_current_ts_extended[i+1] - 1
    a <- cmu[attr[index, i]]
    text(segment_starts_at, x_mat[index, segment_starts_at] + 0.25, round(a, digits=2))
    # text(segment_starts_at, x_mat[index, segment_starts_at] - 0.5, round(a, digits=2))
  }
}

explain_overall <- function(Nattr, Nattr_cp, attr)
{
  cat("Changepoints belief ===>\n")
  Ncp_of_each_ts <- get_number_of_changepoints(cp)
  cp_seatings <- n_customers_at_each_table_cp(Ncp_of_each_ts, Ncp_max)
  names(Ncp_of_each_ts) <- 1:length(Ncp_of_each_ts)
  cat("Number of changepoints for each time series: ", Ncp_of_each_ts, "\n")
  hist(Ncp_of_each_ts)
  cat("Clusters belief ===>\n")
  cat("Number of clusters: ", Nattr, "\n")
  cat("Means of the clusters: ", cmu, "\n")
  cluster_histg <- table(as.vector(attr)[!is.na(attr)])
  names(cluster_histg) <- round(cmu[1:length(cluster_histg)],2)
  barplot(cluster_histg)
}