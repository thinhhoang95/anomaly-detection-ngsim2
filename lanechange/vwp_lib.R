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
    result[i] <- length(cp[i])
  }
  return(result)
}

n_customers_at_each_table_cp <- function(Ncp, max_cp)
{
  result <- rep(0, max_cp)
  t_attr <- table(Ncp)
  for (i in 1:max_cp)
  {
    if (!is.na(t_attr[i]))
    {
      result[i] <- t_attr[i]
    }
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
  cat('Index/Value: ', current_cp_index, current_cp, '\n')
  if (current_cp_index <= max_cp_index)
  {
    current_cp_range <- cp_range(n_length_ts, current_cp_index, current_cp)
    if (length(current_cp_range) > 0)
    {
      if (current_cp <= tail(current_cp_range, n=1))
      {
        
        current_cp <- head(current_cp_range, n=1)
        cat('-> Index/Value: ', current_cp_index, current_cp, '\n')
        
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