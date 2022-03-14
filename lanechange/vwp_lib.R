least_square <- function(x, basis)
{
  # if the basis length is larger than the vector length, just trim off the
  # right part so that it has a length equals to x
  basis <- basis[1:length(x)]
  # the design matrix
  # col_of_ones <- rep(1, length(basis))
  X <- basis
  # the response matrix
  Y <- x
  # get the least square solution
  b <- solve(t(X) %*% X) %*% t(X) %*% Y
  return(list("b" = b, "X" = X))
}

reconstruct <- function(x, b, basis)
{
  basis <- basis[1:length(x)]
  col_of_ones <- rep(1, length(basis))
  X <- cbind(col_of_ones, basis)
  return(X %*% b)
}

cut_x <- function(x_mat, cp, basis, all_pah)
{
  segments <- list()
  as <- list()
  for (ts_index in 1:dim(x_mat)[1])
  {
    segments_of_this_ts <- list()
    as_of_this_ts <- list()
    cp_of_this_ts <- cp[[ts_index]]
    ts_length <- length(x_mat[ts_index,])
    cp_of_this_ts_extended <- c(1, cp_of_this_ts, ts_length + 1)
    for (i in 1:(length(cp_of_this_ts)+1))
    {
      segment_starts_at <- cp_of_this_ts_extended[i]
      segment_ends_at <- cp_of_this_ts_extended[i+1] - 1
      segment <- x_mat[ts_index, segment_starts_at:segment_ends_at]
      segments_of_this_ts <- append(segments_of_this_ts, list(segment))
      a <- all_pah[[ts_index]]$a[segment_starts_at, segment_ends_at]
      as_of_this_ts <- append(as_of_this_ts, list(a))
    }
    segments <- append(segments, list(segments_of_this_ts))
    as <- append(as, list(as_of_this_ts))
  }
  return(list("x_cut" = segments, "a_cut" = as))
}

# To return the segments cut according to cp
cut_x_with_cp <- function(x_cut, cp, basis)
{
  x_cut_result <- list()
  a_cut_result <- list()
  for (ts_index in 1:length(x_cut))
  {
    ts <- x_cut[[ts_index]]
    ts_length = length(ts)
    extended_cp <- append(1, cp[[ts_index]])
    extended_cp <- append(extended_cp, ts_length)
    for (i in 1:(length(extended_cp)-1))
    {
      current_changepoint <- extended_cp[i]
      next_changepoint <- extended_cp[i+1]
      ts_to_append <- ts[current_changepoint:(next_changepoint - 1)]
      segment_a <- as.vector(least_square(ts_to_append, basis)$b)
      if (i==1)
      {
        x_cut_result[ts_index] <- list(ts_to_append) # nothing yet at the row ts_index we are selecting of x_cut_result
        a_cut_result[ts_index] <- list(segment_a)
      }
      else
      {
        x_cut_result[ts_index] <- list(append(list(x_cut_result[[ts_index]]), list(ts_to_append)))
        a_cut_result[ts_index] <- list(append(list(a_cut_result[[ts_index]]), list(segment_a)))
      }
    }
  }
  return(list("x_cut" = x_cut_result, "a_cut" = a_cut_result))
}

n_customers_at_each_table <- function(attr)
{
  result <- rep(0, ncol(attr))
  t_attr <- table(attr)
  t_attr_label <- as.numeric(names(t_attr))
  for (i in 1:length(attr))
  {
    result[t_attr_label[i]] <- t_attr[i]
  }
  return(result)
}

log_p_homogeneity <- function(x, basis, a, vary)
{
  # if the basis length is larger than the vector length, just trim off the
  # right part so that it has a length equals to x
  basis <- basis[1:length(x)]
  tauy <- solve(vary)
  return(-1/2 * tauy * (t(x) %*% x - 2 * t(matrix(x)) %*% matrix(basis) %*% a + t(a) %*% t(basis) %*% basis %*% matrix(a)))
}

get_number_of_changepoints <- function(cp)
{
  result <- rep(0, length(cp))
  for (i in 1:length(cp))
  {
    result[i] <- length(cp[[i]])
  }
  return(result)
}

n_customers_at_each_table_cp <- function(Ncp, max_cp)
{
  result <- rep(0, max_cp)
  t_attr <- table(Ncp)
  t_attr_label <- as.numeric(names(t_attr))
  for (i in 1:length(attr))
  {
    result[t_attr_label[i]] <- t_attr[i]
  }
  return(result)
}

cp_range <- function(n_length_ts, current_cp_index, current_cp)
{
  cp_range <- (current_cp_index + 1) : (n_length_ts - current_cp_index - 1)
  return(cp_range[which(cp_range > current_cp)])
}

# I don't know why but it seems to miss 2 from the n_length_ts. Consider adding 2 when calling this function!
append_cp_recursively <- function(cp_list, vec, n_length_ts, current_cp_index, max_cp_index, current_cp)
{
  # cat('Index/Value: ', current_cp_index, current_cp, '\n')
  if (current_cp_index <= max_cp_index)
  {
    current_cp_range <- cp_range(n_length_ts, current_cp_index, current_cp)
    if (length(current_cp_range) > 0)
    {
      if (current_cp <= tail(current_cp_range, n=1))
      {
        
        current_cp <- head(current_cp_range, n=1)
        # cat('-> Index/Value: ', current_cp_index, current_cp, '\n')
        
        # cat('Index/Value: ', current_cp_index, current_cp, '\n')
        # There are still things to add in this index
        vec[current_cp_index] <- current_cp
        if (current_cp_index == max_cp_index)
        {
          cp_list <- append(cp_list, list(vec))
        }
        cp_list <- append_cp_recursively(cp_list, vec, n_length_ts, current_cp_index + 1, max_cp_index, current_cp)
        # We try to add changepoint current_cp + 1
        cp_list <- append_cp_recursively(cp_list, vec, n_length_ts, current_cp_index, max_cp_index, current_cp)
      }
    }
  }
  return(cp_list)
}

precalculate_a_and_homogeneity <- function(ts, basis, vary)
{
  a_array <- array(, dim=c(length(ts), length(ts)))
  p_array <- array(, dim=c(length(ts), length(ts)))
  for (i in 1:(length(ts)-1))
  {
    for (j in i:length(ts))
    {
      # Calculate the slope
      a <- as.vector(least_square(ts[i:j], basis)$b)
      a_array[i,j] <- a 
      # Calculate the homogeneity log likelihood
      p <- log_p_homogeneity(ts[i:j], basis, a, vary)
      p_array[i,j] <- p
    }
  }
  return(list("a" = a_array, "p" = p_array))
}

find_most_likely_partitioning_of_a_time_series <- function(ts_length, Ncp, cmu, ctau, Nattr, precalculated_ah)
{
  # ts is the vector of time series
  # Ncp is the number of changepoints one allowed to have
  # cmu is the cluster mean
  # cvar is considered to be a constant for all clusters
  
  # Get all possible changepoint placements
  cp_placements <- append_cp_recursively(list(), rep(0, Ncp), ts_length + 2, 1, Ncp, 0) # the +2 is just a workaround!!! Remember to +2 onto every ts length
  p_cp_placement <- rep(0, length(cp_placements))
  optimal_cluster_attr <- list()
  for (l in 1:length(cp_placements))
  {
    cp_placement <- cp_placements[[l]]
    cp_placement_extended <- c(1, cp_placement, ts_length + 1)
    for (i in 1:(length(cp_placement)+1))
    {
      segment_starts_at <- cp_placement_extended[i]
      segment_ends_at <- cp_placement_extended[i+1] - 1
      # We get a and homogeneity loglikelihood for the corresponding segment
      a <- precalculated_ah$a[segment_starts_at, segment_ends_at]
      # Get the most likely cluster that a belongs to 
      p_a_in_cluster <- rep(0, Nattr)
      for (cluster_id in 1:Nattr)
      {
        p_a_in_cluster[cluster_id] <- -1/2 * ctau * (a - cmu[cluster_id])^2
      }
      lp_homo <- precalculated_ah$p[segment_starts_at, segment_ends_at]
      p_cp_placement[l] <- p_cp_placement[l] + lp_homo + max(p_a_in_cluster)
      if (i==1)
      {
        best_cluster_for_this_segment <- which.max(p_a_in_cluster)
      } else 
      {
        best_cluster_for_this_segment <- c(best_cluster_for_this_segment, which.max(p_a_in_cluster))
      }
    }
    optimal_cluster_attr <- append(optimal_cluster_attr, list(best_cluster_for_this_segment))
  }
  index_of_best_changepoints <- which.max(p_cp_placement)
  return(list("cp" = cp_placements[[index_of_best_changepoints]], "lp" = p_cp_placement[index_of_best_changepoints],
              "attr" = optimal_cluster_attr[[index_of_best_changepoints]]))
}