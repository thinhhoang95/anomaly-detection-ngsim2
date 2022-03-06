least_square <- function(x, basis)
{
  # if the basis length is larger than the vector length, just trim off the
  # right part so that it has a length equals to x
  basis <- basis[1:length(x)]
  # the design matrix
  col_of_ones <- rep(1, length(basis))
  X <- cbind(col_of_ones, basis)
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
cut_x_with_cp <- function(x_cut, cp)
{
  x_cut_result <- list()
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
      if (i==1)
      {
        x_cut_result[ts_index] <- list(ts_to_append) # nothing yet at the row ts_index we are selecting of x_cut_result
      }
      else
      {
        x_cut_result[ts_index] <- list(append(list(x_cut_result[[ts_index]]), list(ts_to_append)))
      }
    }
  }
  return(x_cut_result)
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