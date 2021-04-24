#' Function that computes the Maximum Likelihood PCA with equal row correlated errors.
#'
#' @author Tobias K. Karakach, Federico Taverna
#'
#' @import data.table
#'
#' @param expressionMatrix The expression matrix organized as (*n* X *m*), *n* is the number of genes and *m* is the number of samples.
#' You should use the expressionMatrix output of prepareData().
#'
#' @param errorCovarianceMatrix The inverse of the measurement errors, whose dimensions are (*n* X *n*).
#' You should use the errorCovarianceMatrix output of prepareData().
#'
#' @param numberOfComponents The number of principal components estimated to be required to reconstruct the expression matrix without loss of generality.
#'
#' @param maxIterations The maximum number of iterations to peform during the computation of the Maximum Likelihood. Once the max iteration is reached,
#' the algorithms stops.
#'
#' @param tolerance The termination tolerance for the computation of the Maximum Likelihood. Once the improvements made by Maximum Likelihood reach the
#' tolerance, the algorithm stops.
#'
#' @param verbose Whether to display the information about the computation or not.
#'
#' @return A list consisting of U (principal components), S(eigenvalues) and V (loadings).
#' The estimated matrix (estimatedMatrix) is given as t(U %*% S %*% t(V)).
#'
#' @export
mlProjection <- function(expressionMatrix, errorCovarianceMatrix, numberOfComponents = 15, maxIterations = 2000, tolerance = 1e-10, verbose = TRUE)
{
  # Start timer
  start <- Sys.time()

  # Sanity check for the user input
  if (!any("matrix" %in% class(expressionMatrix))) stop ("The expression matrix must be in the matrix format.")
  if (!any("matrix" %in% class(errorCovarianceMatrix))) stop ("The inverse error covariance matrix must be in the matrix format.")
  if (!"numeric" %in% class(numberOfComponents)) stop("The number of components must be a number.")
  if (!"logical" %in% class(verbose)) stop("The verbose variable must be logical (TRUE or FALSE)")

  # Assign the arguments to new variable names
  X <- expressionMatrix
  Xcov <- errorCovarianceMatrix
  ncomp <- numberOfComponents
  rm(expressionMatrix, errorCovarianceMatrix, numberOfComponents)

  # Transpose the expression matrix
  n <- dim(X)[1]
  m <- dim(X)[2]
  X <- t(X)

  # Set the maximum likelihood parameters
  pc <- ncomp
  flg <- 1
  maxiter <- maxIterations
  lamda <- tolerance
  iter <- 0

  # Compute the singular-value decomposition
  tmp <- svd(X)
  u <- tmp$u
  v <- tmp$v
  s <- diag(tmp$d)

  # Select v with the right number of principal components
  V <- tmp$v[, 1:pc]

  # Maximum likelihood loop
  while (flg == 1 && iter < maxiter)
  {
    # Increment iteration count
    iter <- iter + 1

    # Comment TBD
    Xhat <- X %*% ((Xcov)) %*% V %*% solve(t(V) %*% (Xcov) %*% V) %*% t(V) # Estimated matrix
    S1 <- sum(diag(((X - Xhat) %*% (Xcov) %*% t((X - Xhat)))))

    # Comment TBD
    tmp <- svd(Xhat)
    Xhat <- tmp$u[, 1:pc] %*% t(tmp$u[,1:pc]) %*% X  # Estimated matrix
    S2 <- sum(diag(((X - Xhat) %*% (Xcov) %*% t((X-Xhat)))))
    tmp <- svd(Xhat)
    U <- tmp$u[, 1:pc]
    S <- diag(tmp$d[1:pc])
    V <- tmp$v[, 1:pc]

    # If the threshold is reached, get out of the loop
    sobj <- abs((S1-S2)/S2)
    if(sobj < lamda) flg <- 0

    # Print the iteration information
    if (verbose == TRUE) message(paste0("Iteration number: ", iter, "; Converging to: ", lamda, "; Currently at: ", sobj))
  }

  # Assign the correct names
  rownames(U) <- rownames(X)
  colnames(U) <- paste0("Comp", 1:ncomp)
  rownames(V) <- colnames(X)
  colnames(V) <- paste0("Comp", 1:ncomp)
  rownames(S) <- paste0("Comp", 1:ncomp)
  colnames(S) <- paste0("Comp", 1:ncomp)

  # Compute the estimated matrix
  estimatedMatrix <- U %*% S %*% t(V)
  estimatedMatrix <- t(estimatedMatrix)

  # End timer and display elapsed time
  end <- Sys.time()
  minutes <- round(difftime(end, start, units = "mins"), digits = 2)
  if (verbose == TRUE) message(paste0("ML projection computed in ", minutes, " minutes."))

  # Create the output list
  outputList <- list("U" = U, "S" = S, "V" = V, "estimatedMatrix" = estimatedMatrix)

  return(outputList)
}
