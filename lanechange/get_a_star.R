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