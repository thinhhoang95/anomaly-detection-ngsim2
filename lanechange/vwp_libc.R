# This version will try to ensure the continuity of the reconstructed curve
# Consider merging the clusters

least_square <- function(x_original, basis)
{
  # if the basis length is larger than the vector length, just trim off the
  # right part so that it has a length equals to x
  basis <- basis[1:length(x_original)]
  # the design matrix
  # col_of_ones <- rep(1, length(basis))
  X <- basis
  # the response matrix
  # make sure the ts starts from zero to get the slope
  x <- x_original - x_original[1]
  Y <- x
  # get the least square solution
  b <- solve(t(X) %*% X) %*% t(X) %*% Y
  return(list("b" = b, "X" = X))
}

reconstruct <- function(x_length, b, basis)
{
  basis <- basis[1:x_length]
  #col_of_ones <- rep(1, length(basis))
  #X <- cbind(col_of_ones, basis)
  return(b %*% basis)
}

# This is the legacy version that uses the homogeneity likelihood instead of the joint likelihood
# see the cut_x function below
cut_x_legacy_pah <- function(x_mat, cp, basis, all_pah)
{
  segments <- list()
  as <- list()
  for (ts_index in 1:dim(x_mat)[1])
  {
    segments_of_this_ts <- list()
    as_of_this_ts <- list()
    cp_of_this_ts <- cp[[ts_index]]
    ts_length <- length(x_mat[ts_index,])
    cp_of_this_ts_extended <- c(1, cp_of_this_ts, ts_length)
    for (i in 1:(length(cp_of_this_ts)+1))
    {
      segment_starts_at <- cp_of_this_ts_extended[i]
      segment_ends_at <- cp_of_this_ts_extended[i+1]
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

cut_x <- function(x_mat, cp, basis, pa_homo)
{
  segments <- list()
  as <- list()
  for (ts_index in 1:dim(x_mat)[1])
  {
    segments_of_this_ts <- list()
    as_of_this_ts <- list()
    cp_of_this_ts <- cp[[ts_index]]
    ts_length <- length(x_mat[ts_index,])
    cp_of_this_ts_extended <- c(1, cp_of_this_ts, ts_length)
    for (i in 1:(length(cp_of_this_ts)+1))
    {
      segment_starts_at <- cp_of_this_ts_extended[i]
      segment_ends_at <- cp_of_this_ts_extended[i+1]
      segment <- x_mat[ts_index, segment_starts_at:segment_ends_at]
      segments_of_this_ts <- append(segments_of_this_ts, list(segment))
      a <- pa_homo$a[ts_index, segment_starts_at, segment_ends_at]
      as_of_this_ts <- append(as_of_this_ts, list(a))
    }
    segments <- append(segments, list(segments_of_this_ts))
    as <- append(as, list(as_of_this_ts))
  }
  return(list("x_cut" = segments, "a_cut" = as))
}

# To return the segments cut according to cp
# Attention: this function is deprecated at the moment, replaced by the
# cut_x function above
# cut_x_with_cp <- function(x_cut, cp, basis)
# {
#   x_cut_result <- list()
#   a_cut_result <- list()
#   for (ts_index in 1:length(x_cut))
#   {
#     ts <- x_cut[[ts_index]]
#     ts_length = length(ts)
#     extended_cp <- append(1, cp[[ts_index]])
#     extended_cp <- append(extended_cp, ts_length)
#     for (i in 1:(length(extended_cp)-1))
#     {
#       current_changepoint <- extended_cp[i]
#       next_changepoint <- extended_cp[i+1]
#       ts_to_append <- ts[current_changepoint:(next_changepoint - 1)]
#       segment_a <- as.vector(least_square(ts_to_append, basis)$b)
#       if (i==1)
#       {
#         x_cut_result[ts_index] <- list(ts_to_append) # nothing yet at the row ts_index we are selecting of x_cut_result
#         a_cut_result[ts_index] <- list(segment_a)
#       }
#       else
#       {
#         x_cut_result[ts_index] <- list(append(list(x_cut_result[[ts_index]]), list(ts_to_append)))
#         a_cut_result[ts_index] <- list(append(list(a_cut_result[[ts_index]]), list(segment_a)))
#       }
#     }
#   }
#   return(list("x_cut" = x_cut_result, "a_cut" = a_cut_result))
# }

n_customers_at_each_table <- function(attr, Nattr)
{
  t_attr <- table(attr)
  t_attr_label <- as.numeric(names(t_attr))
  N_clusters_to_store = max(c(Nattr, max(t_attr_label)))
  result <- rep(0, N_clusters_to_store)
  for (i in 1:N_clusters_to_store)
  {
    result[t_attr_label[i]] <- t_attr[i]
  }
  return(result)
}

log_p_homogeneity <- function(x_original, basis, a, vary)
{
  # if the basis length is larger than the vector length, just trim off the
  # right part so that it has a length equals to x
  basis <- basis[1:length(x_original)]
  x <- x_original - x_original[1]
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

precalculate_pa_homo <- function(x_mat, basis, vary)
{
  result_p <- array(NA, dim=c(nrow(x_mat), ncol(x_mat), ncol(x_mat)))
  result_a <- array(NA, dim=c(nrow(x_mat), ncol(x_mat), ncol(x_mat)))
  for (ts_index in 1:nrow(x_mat))
  {
    ts <- x_mat[ts_index,]
    # so the time series ts would start from zero, not from any arbitrary value
    a_array <- array(, dim=c(length(ts), length(ts)))
    p_array <- array(, dim=c(length(ts), length(ts)))
    for (i in 1:(length(ts)-1))
    {
      for (j in i:length(ts))
      {
        # Calculate the slope
        if (j==i)
        {
          # In this case, the ts contains only {0} so the slope is undetermined
          # We can choose an arbitrary slope, let's say 0
          a <- 0
          p <- -Inf # we do not want the algorithm to cheat by using an infinite amount of changepoints to get zero probability
        } else if (j==i+1)
        {
          # In this case, the slope can be estimated but nevertheless, the log proba is still zero (i.e. proba=1)
          a <- as.vector(least_square(ts[i:j], basis)$b)
          p <- -Inf # we try to avoid this also
        } else {
          a <- as.vector(least_square(ts[i:j], basis)$b)
          # Calculate the homogeneity log likelihood
          p <- log_p_homogeneity(ts[i:j], basis, a, vary)
        }
        a_array[i,j] <- a 
        p_array[i,j] <- p
      }
    }
    result_p[ts_index,,] <- p_array 
    result_a[ts_index,,] <- a_array
  } # rows of x_mat
  return(list("a" = result_a, "p" = result_p))
}

# This will find the likelihood of every substring of the timeseries, marginalized by the clusters
# this is different from the precalculated_a_and_homogeneity, which only deals with the homogeneity of the time series
# not how typicality will the ts fit into the cluster
precalculate_a_and_likelihood <- function(ts_index, ts_length, pa_homo, basis, vary, Nattr, cmu, ctau, segment_seatings, alpha)
{
  # so the time series ts would start from zero, not from any arbitrary value
  a_array <- array(, dim=c(ts_length, ts_length))
  p_array <- array(, dim=c(ts_length, ts_length))
  for (i in 1:(ts_length-1))
  {
    for (j in i:ts_length)
    {
      # Calculate the slope
      if (j==i)
      {
        # In this case, the ts contains only {0} so the slope is undetermined
        # We can choose an arbitrary slope, let's say 0
        a <- 0
        p <- -Inf # we do not want the algorithm to cheat by using an infinite amount of changepoints to get zero probability
      } else if (j==i+1)
      {
        # In this case, the slope can be estimated but nevertheless, the log proba is still zero (i.e. proba=1)
        # a <- as.vector(least_square(ts[i:j], basis)$b)
        a <- pa_homo$a[ts_index, i, j]
        p <- -Inf # we try to avoid this also
      } else {
        a <- pa_homo$a[ts_index, i, j]
        # a <- as.vector(least_square(ts[i:j], basis)$b)
        # Calculate the expected typicality probability
        p_a_in_cluster <- rep(0, Nattr)
        p_typicality_in_cluster <- rep(0, Nattr)
        for (cluster_id in 1:Nattr)
        {
          # p_a_in_cluster[cluster_id] <- -1/2 * ctau * (a - cmu[cluster_id])^2
          p_a_in_cluster[cluster_id] <- dnorm(a, mean = cmu[cluster_id], sd=sqrt(solve(ctau)))
          nk <- segment_seatings[cluster_id]
          total_seatings <- sum(segment_seatings)
          p_popularity_of_the_current_cluster <- nk/(total_seatings - 1 + alpha) # the mixture proportion
          p_typicality_in_cluster[cluster_id] <- p_popularity_of_the_current_cluster * p_a_in_cluster[cluster_id]
        }
        p_expected_typicality_among_all_clusters <- sum(p_typicality_in_cluster)
        
        # Calculate the homogeneity log likelihood
        p_homogeneity <- pa_homo$p[ts_index, i, j]
        # p_homogeneity <- log_p_homogeneity(ts[i:j], basis, a, vary)
        # The proba of the time series ij, expected over the current clusters:
        p <- log(p_expected_typicality_among_all_clusters) + p_homogeneity
      }
      a_array[i,j] <- a 
      p_array[i,j] <- p
    }
  }
  return(list("a" = a_array, "p" = p_array))
}

precalculate_a_and_homogeneity <- function(ts, basis, vary)
{
  # so the time series ts would start from zero, not from any arbitrary value
  a_array <- array(, dim=c(length(ts), length(ts)))
  p_array <- array(, dim=c(length(ts), length(ts)))
  for (i in 1:(length(ts)-1))
  {
    for (j in i:length(ts))
    {
      # Calculate the slope
      if (j==i)
      {
        # In this case, the ts contains only {0} so the slope is undetermined
        # We can choose an arbitrary slope, let's say 0
        a <- 0
        p <- -Inf # we do not want the algorithm to cheat by using an infinite amount of changepoints to get zero probability
      } else if (j==i+1)
      {
        # In this case, the slope can be estimated but nevertheless, the log proba is still zero (i.e. proba=1)
        a <- as.vector(least_square(ts[i:j], basis)$b)
        p <- -Inf # we try to avoid this also
      } else {
        a <- as.vector(least_square(ts[i:j], basis)$b)
        # Calculate the homogeneity log likelihood
        p <- log_p_homogeneity(ts[i:j], basis, a, vary)
      }
      a_array[i,j] <- a 
      p_array[i,j] <- p
    }
  }
  return(list("a" = a_array, "p" = p_array))
}

### THIS IS THE LEGACY METHOD TO FIND THE CHANGEPOINTS THAT CAN BE EXTREMELY SLOW WITH LARGE NUMBER OF CHANGEPOINTS
### USE THE NEW METHOD: DYNAMIC PROGRAMMING INSTEAD
find_most_likely_partitioning_of_a_time_series <- function(ts_length, Ncp, cmu, ctau, Nattr, precalculated_ah)
{
  # ts is the vector of time series
  # Ncp is the number of changepoints one allowed to have
  # cmu is the cluster mean
  # cvar is considered to be a constant for all clusters
  
  # Get all possible changepoint placements
  cp_placements <- append_cp_recursively(list(), rep(0, Ncp), ts_length + 2, 1, Ncp, 0) # the +2 is just a workaround!!! Remember to +2 onto every ts length
  # cat("Placements obtained! \n")
  p_cp_placement <- rep(0, length(cp_placements))
  optimal_cluster_attr <- list()
  segment_slopes_for_each_cp_placement <- list()
  for (l in 1:length(cp_placements)) # loop through the CP placements
  {
    cp_placement <- cp_placements[[l]]
    cp_placement_extended <- c(1, cp_placement, ts_length)
    for (i in 1:(length(cp_placement)+1)) # loop through the segments
    {
      segment_starts_at <- cp_placement_extended[i]
      segment_ends_at <- cp_placement_extended[i+1]
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
        a_for_this_cp_placement <- a
      } else 
      {
        best_cluster_for_this_segment <- c(best_cluster_for_this_segment, which.max(p_a_in_cluster))
        a_for_this_cp_placement <- c(a_for_this_cp_placement, a)
      }
    }
    optimal_cluster_attr <- append(optimal_cluster_attr, list(best_cluster_for_this_segment))
    segment_slopes_for_each_cp_placement <- append(segment_slopes_for_each_cp_placement, list(a_for_this_cp_placement))
  }
  index_of_best_changepoints <- which.max(p_cp_placement)
  return(list("cp" = cp_placements[[index_of_best_changepoints]], "lp" = p_cp_placement[index_of_best_changepoints],
              "attr" = optimal_cluster_attr[[index_of_best_changepoints]], "a" = segment_slopes_for_each_cp_placement[[index_of_best_changepoints]]))
}

dynamic_program_cp <- function(sol, Ncp, ts, ts_index, pa_homo, k, basis, vary, Nattr, cmu, ctau, segment_seatings, alpha)
{
  # Use dynamic programming to fill in the solution of finding changepoints, shown in the sol matrix
  # default_sol <- list("cp" = array(,dim=c(10,20,10)), "p" = array(,dim=c(20,10)))
  # 20: length of ts, 10: maximum number of changepoints
  if (k>Ncp)
  {
    return(sol)
  }
  
  ts_length <- length(ts)
  # ts_pah <- precalculate_a_and_homogeneity(ts, basis, vary) # only has p_homogeneity,
  ts_pah <- precalculate_a_and_likelihood(ts_index, ts_length, pa_homo, basis, vary, Nattr, cmu, ctau, segment_seatings, alpha)
  for (v in 1:ts_length)
  {
    sub_ts <- ts[v:ts_length]
    if (v+1 > ts_length - k - 1)
    {
      # The sub_ts will be too short to contain a meaningful changepoint, will do nothing
    } else {
      range_of_cp <- (v+1):(ts_length-k-1)
      logp_of_subts_with_each_cp_placement <- rep(0, length(range_of_cp))
      a_of_the_first_segment <- rep(0, length(range_of_cp))
      a_of_the_second_segment <- list()
      for (proposed_cp in range_of_cp)
      {
        first_segment_starts_at <- v
        first_segment_ends_at <- proposed_cp 
        second_segment_starts_at <- proposed_cp 
        second_segment_ends_at <- ts_length
        
        if (k==1)
        {
          # This is the only changepoint left in the time series
          logp_second_segment <- ts_pah$p[second_segment_starts_at, second_segment_ends_at]
          a_of_the_second_segment <- append(a_of_the_second_segment, list(ts_pah$a[second_segment_starts_at, second_segment_ends_at]))
        } else { # k>1
          logp_second_segment <- sol$p[second_segment_starts_at,k-1]
          a_of_the_second_segment <- append(a_of_the_second_segment, list(sol$a[1:k, second_segment_starts_at, k-1])) # copies the a as well
        }
        
        logp_of_subts_with_each_cp_placement[proposed_cp - v] <- ts_pah$p[first_segment_starts_at, first_segment_ends_at] + logp_second_segment
        a_of_the_first_segment[proposed_cp - v] <- ts_pah$a[first_segment_starts_at, first_segment_ends_at]
      }
      index_of_cp_that_yields_max_logp <- which.max(logp_of_subts_with_each_cp_placement)
      cp_that_yields_max_logp <- range_of_cp[index_of_cp_that_yields_max_logp]
      max_logp <- logp_of_subts_with_each_cp_placement[index_of_cp_that_yields_max_logp]
      best_a_of_the_first_segment <- a_of_the_first_segment[index_of_cp_that_yields_max_logp]
      best_a_of_the_second_segment <- unlist(a_of_the_second_segment[[index_of_cp_that_yields_max_logp]])
      sol$p[v,k] <- max_logp
      sol$a[1,v,k] <- best_a_of_the_first_segment
      sol$a[2:(k+1),v,k] <- best_a_of_the_second_segment
      sol$cp[1,v,k] <- cp_that_yields_max_logp # the first changepoint
      if (k>1)
      {
        sol$cp[2:k,v,k] <- sol$cp[1:k-1, cp_that_yields_max_logp, k-1] # other changepoints are copied
        # from the solution to the time series assumed to start from the just found changepoint
      }
    } # check if the sub_ts is too short
    
  } # loop through v (the substring beginning index)
  if (k < Ncp)
  {
    # Recursive / dynamic programming
    sol <- dynamic_program_cp(sol, Ncp, ts, ts_index, pa_homo, k+1, basis, vary, Nattr, cmu, ctau, segment_seatings, alpha)
  }
  return(sol)
}

merge_segments <- function(cp, attr)
{
  # This function will merge adjacent segments with similar cluster attribution
  # and will return the new number of changepoints
  for (row in 1:length(cp))
  {
    # For each time series (row), we run through the CP associated
    attr_row <- attr[row,]
    attr_row_without_na <- attr[row,][!is.na(attr[row,])]
    n_segments <- length(attr_row_without_na)
    if (n_segments > 1)
    {
      for (col in 1:(n_segments - 1))
      {
        if (attr_row_without_na[col] == attr_row_without_na[col+1])
        {
          # We find an attr like 3,3 so it is wise to remove the changepoint associated
          index_of_cp_to_remove <- col # the index of the changepoint to remove would be the same as col
          cp_at_row <- cp[[row]]
          new_cp_at_row <- cp_at_row[-index_of_cp_to_remove]
          if (length(new_cp_at_row) <= 1)
          {
            # Do nothing, since we don't want the time series to have zero changepoints
          } else {
            # Only merge segments if we have more than 1 changepoint left
            cp[row] <- list(new_cp_at_row)
          }
          # We also update the attr(ibution)
          # The idea is to delete the entry at col+1 and append a NA at the end
          attr_row <- c(attr_row[-(col+1)], NA)
          attr[row,] <- attr_row
          recufunc <- merge_segments(cp, attr)
          # The return command will essentially break all loops
          # The reason why we did this was because cp and attr were changed, so looping further on them
          # would cause problems with indices
          return(list("cp" = recufunc$cp, "attr" = recufunc$attr))
        }
      }
    }
  }
  return(list("cp" = cp, "attr" = attr))
}

get_Nattr_cp_from_attr <- function(attr)
{
  segments_in_ts <- rep(0, nrow(attr))
  for (row in 1:nrow(attr))
  {
    attr_row <- attr[row,]
    attr_row <- attr_row[!is.na(attr_row)] # remove all NAs
    segments_in_ts[row] <- length(attr_row)
  }
  return(max(segments_in_ts))
}