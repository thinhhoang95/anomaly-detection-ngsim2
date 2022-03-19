visualize_belief <- function(x_mat, cp, cmu, attr, index)
{
  ts_length <- ncol(x_mat)
  cp_of_current_ts <- cp[[index]]
  cp_of_current_ts_extended <- c(1, cp_of_current_ts, ts_length + 1)
  x_recons <- rep(0, ts_length)
  for (i in 1:(length(cp_of_current_ts)+1))
  {
    segment_starts_at <- cp_of_current_ts_extended[i]
    segment_ends_at <- cp_of_current_ts_extended[i+1] - 1
    cat("\n Segment starts at", segment_starts_at)
    cat("\n Segment ends at", segment_ends_at)
    # Reconstruct the segment
    a <- cmu[attr[index, i]]
    for (t in segment_starts_at:segment_ends_at)
    {
      if (t==1)
      {
        x_recons[t] <- 0
      } else {
        x_recons[t] <- x_recons[t-1] + a
      }
    }
  }
  
  plot(x_recons, col='blue')
  points(x_mat[index,])
}