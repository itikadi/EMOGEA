#' Function to perform multivariate curve resolution by alternating total least squares.
#'
#' @author Tobias K. Karakach, Federico Taverna
#'
#' @import data.table
#' @import Matrix
#'
#' @param expressionMatrix The expression matrix organized as (*n* X *m*), *n* is the number of genes and *m* is the number of samples.
#' You should use the expressionMatrix output of prepareData().
#'
#' @param residualMatrix The residual matrix, whose dimensions are (*n* X *m*). You should use the residualMatrix output of prepareData().
#' If the residual matrix is specified, the function uses the error weighted multivariate curve resolution algorithm.
#' If the residual matrix is NULL, the function uses the non-weighted multivariate curve resolution algorithm.
#'
#' @param numberOfComponents The number of number of components (profiles) to be extracted.
#'
#' @param initAlgorithm Which algorithm to use to initialize the P/C values. The choice is between "simplisma" (recommended) and "random".
#'
#' @param maxIterations The maximum number of iteration for the maximum likelihood algorithm.
#'
#' @param tolerance The convergence tolerance. Once this is reached, the algorithm stops.
#'
#' @param randomSeed The random seed number, only used if initAlgorithm is "random".
#'
#' @param computeGeneProfiles Whether to compute the gene profiler or not.
#'
#' @param verbose Whether to display the information about the computation or not.
#'
#' @return The output of is a list containing: Popt, the final estimate of the profile (spectral) matrix, where
#' vectors (rows) are constrained to be non-negative and unit length. Copt, the final estimate of contribution matrix,
#' constrained to be non-negative. estimatedMatrix, the final estimated matrix. geneProfiles, the gene profiles (only if computeGeneProfiles = TRUE).
#'
#' @export
multivariateCurveResolution <- function(expressionMatrix, residualMatrix = NULL, numberOfComponents = 3,
                                        initAlgorithm = "simplisma", maxIterations = 2000, tolerance = 0.0009,
                                        randomSeed = 091120, computeGeneProfiles = TRUE, verbose = TRUE)
{
  # Start timer
  start <- Sys.time()

  # Assign the arguments to new variable names
  X <- expressionMatrix
  Xres <- residualMatrix
  ncomp <- numberOfComponents
  maxiter <- maxIterations
  rm(expressionMatrix, residualMatrix, numberOfComponents, maxIterations)

  # Choose algorithm based on whether the residual matrix is given as input
  if (!is.null(Xres)) alg <- "weighted"
  else alg <- "nonweighted"

  # Get the expression matrix dimensions (n = genes,  m = samples)
  n <- dim(X)[1]
  m <- dim(X)[2]

  # ----------------------------------------------------
  # Initialization section
  # ----------------------------------------------------
  tol <- tolerance # Set tolerance for convergence
  flg <- 0 # This is set to 1 to stop the iterative procedure

  # Transpose the expression matrix and residual matrix
  Xoriginal <- X
  X <- t(X)
  if (!is.null(Xres)) Xres <- t(Xres)

  # Initialize P
  if(initAlgorithm == "simplisma")
  {
    simplismaOutput <- simplisma(expressionMatrix = Xoriginal, numberOfComponents = ncomp)
    Prof <- simplismaOutput$Pinit
  }
  else if (initAlgorithm == "random")
  {
    set.seed(randomSeed)       # Set the random seed for the sampling
    itmp <- sample(m)          # Measurement vectors
    Prof <- X[itmp[1:ncomp], ] # Starting points for the iteration from random samples
    Prof[(Prof < 0)] <- 0      # Force any negative values in initial to zero
  }

  # Define the euclidean distance function
  euc.dist <- function(x) { sqrt(sum((x - mean(x)) ^ 2)) }

  # Define the root mean square deviation
  rmse <- function (actual, predicted) { return(sqrt((mean(actual - predicted) ^ 2))) }

  # Normalize profiles to unit length
  for (i in 1:dim(Prof)[1])
  {
    Prof[i,] <- Prof[i, ] / euc.dist(Prof[i, ])
  }

  # ----------------------------------------------------
  # Proceed with alternating TLS algorithm
  # ----------------------------------------------------
  P <- Prof # Starting estimates

  if (alg == "weighted")
  {
    Xwt <- 1 / (Xres ^ 2)   # Weighting factors are reciprocals of variances
    icnt <- 0               # Initialize iteration counter
    rmsdif <- 1

    # Loop to carry out alternating orthogonal LS iterations
    while (rmsdif > tol & icnt < maxiter)
    {
      Pold <- P         # Save old profile matrix
      icnt <- icnt + 1  # Increment iteration count

      # Display iteration count (if necessary)
      if (verbose == TRUE) message(paste0('Alternating Total Least Squares iteration ', icnt))

      # Initialize matrix of calculated data
      Xcalc <- matrix(0, nrow = m, ncol = n, dimnames = list(rownames(X), colnames(X)))

      # Loop to do a maximum likelihood projection of each row onto the subspace of P
      for (i in 1:m)
      {
        values <- Xwt[i, ]
        QQ <- Matrix::Diagonal(n = length(values), x = values)
        Xtrunc <- P %*% QQ %*% t(P)
        Xstrt <- X[i, ] %*% QQ %*% t(P)
        rm(QQ)

        # Maximum likelihood projection of row i of X onto estimated P
        Xcalc[i, ] <-  as.matrix(Xstrt %*% solve(Xtrunc) %*% P)
        rm(Xstrt, Xtrunc)
      }

      C <- Xcalc %*% t(P) %*% solve(P %*% t(P)) # Solve for C using projected data.
      C[C < 0] <- 0 # Force values less than zero to 0
      Cold <- C

      # Loop to normalize vectors in the contribution matrix
      for (i in 1:ncomp)
      {
        C[, i] <- C[, i]/ euc.dist(C[, i])
      }

      # Re-initialize calculated X for projection into alternate space
      Xcalc <- matrix(0, nrow = m, ncol = n, dimnames = list(rownames(X), colnames(X)))

      # Loop to perform a maximum likelihood
      for (i in 1:n)
      {
        # Projection of each column onto subspace of C.
        # values <- Xwt[, i]
        # QQ <- Matrix::Diagonal(n = length(values), x = values)
        QQ <- diag(Xwt[, i]) # Do not use sparse matrix, as it slows it down too much in this case

        # QQ is the diagonal matrix representing the inverse of the column error covariance matrix
        # Maximum likelihood projection of column i of X onto estimated C.
        Xcalc[, i] <- as.matrix(C %*% solve(t(C) %*% QQ %*% C) %*% (t(C) %*% (QQ %*% X[, i])))
      }

      # Solve for P using projected data.
      P <- solve(t(C) %*% C) %*% (t(C)) %*% Xcalc
      P[P < 0] <- 0

      # Loop to normalize vectors in the profile matrix.
      for (i in 1:ncomp)
      {
        P[i, ] <- P[i, ]/ euc.dist(P[i, ])
      }

      #-------------------------------------------------------------------
      # One iteration complete - now check for convergence. This version
      # uses an absolute tolerance on the root-mean-squared differences
      # between successive normalized profiles
      #-------------------------------------------------------------------
      rms <- rmse(P, Pold)
      if (rms > rmsdif)
      {
        flg <- 1
        if (verbose == TRUE) message("Warning! Possible local minimum. Consider different initial values.")
        warning("Possible local minimum. Consider different initial values.")
      }

      rmsdif <- rms
      if (verbose == TRUE) message(paste0("Converging to: ", tol, "; Currently at: ", rms))

      # Calculate rms differences
      # Check to see if smallest difference is below tolerance
      if (rmsdif < tol)
      {
        flg <- 1
        if (verbose == TRUE) message("Converged. Calculating C with ML projections of X onto P.")
      }

      # Check to see if max iterations exceeded
      if (icnt == maxiter)
      {
        flg <- 1
        warning("Maximum iterations exceeded - stopped with no convergence. Calculating C with ML projections of X onto P.")
      }
    } # End of while

    # Initialize calculated X values
    Xcalc <- matrix(0, nrow = m, ncol = n, dimnames = list(rownames(X), colnames(X)))

    # Loop to do a maximum likelihood projection of each row onto the subspace of P
    for (i in 1:dim(Xcalc)[1])
    {
      values <- Xwt[i, ]
      QQ <- Matrix::Diagonal(n = length(values), x = values)
      Xtrunc <- P %*% QQ %*% t(P)
      Xstrt <- X[i, ] %*% QQ %*% t(P)
      rm(QQ)

      # Maximum likelihood projection of row i of X onto estimated P.
      Xcalc[i, ] <- as.matrix(Xstrt %*% solve(Xtrunc) %*% P)
    }

    C <- Xcalc %*% t(P) %*% solve(P %*% t(P)) # Solve for C using projected data
    C[C < 0] <-0                              # Force nagative values to zero

    # Transpose Xcalc
    Xcalc <- t(Xcalc)

  } # End of weighted

  if (alg == "nonweighted")
  {
    icnt <- 0
    rmsdif <- 1

    while (rmsdif > tol & icnt < maxiter)
    {
      Pold <- P         # Save old profile matrix
      icnt <- icnt + 1  # Increment iteration count

      # Display iteration count (if necessary)
      if (verbose == TRUE) message(paste0("Alternating Least Squares iteration ", icnt))

      C <- X %*% t(P) %*% solve(P %*% t(P)) # Solve for C using projected data
      C[C < 0] <- 0                         # Force values less than zero to 0
      Cold <- C

      # Loop to normalize vectors in the contribution matrix
      for (i in 1:ncomp)
      {

        C[, i] <- C[, i] / euc.dist(C[, i])
      }

      # Solve for P using projected data
      P <- solve(t(C) %*% C) %*% (t(C)) %*% X
      P[P < 0] <- 0

      # Loop to normalize vectors in the profile matrix
      for (i in 1:ncomp)
      {

        P[i, ] <- P[i, ] / euc.dist(P[i, ])
      }

      #----------------------------------------------------------------
      # One iteration complete - now check for convergence. This version
      # uses an absolute tolerance on the root-mean-squared differences
      # between successive normalized profiles.
      #-------------------------------------------------------------------
      rmsdif <- rmse(P, Pold)
      if (verbose == TRUE) message(paste0("Converging to: ", tol, "; Currently at: ", rmsdif))

      # Calculate rms differences
      # Check to see if smallest difference is below tolerance
      if (rmsdif < tol | icnt == maxiter)
      {
        flg <- 1
        if (verbose == TRUE) message('Converged.')
      }

      # Check to see if max iterations exceeded
      if (icnt == maxiter)
      {
        flg <- 1
        warning('Maximum iterations exceeded - stopped with no convergence.')
      }
    } # End of while

    # Calculate estimated matrix
    Xcalc <- t(C %*% P)
    rownames(Xcalc) <- rownames(Xoriginal)
    colnames(Xcalc) <- colnames(Xoriginal)

  } # End of nonweighted

  geneProfiles <- NULL
  if (computeGeneProfiles == TRUE)
  {
    if (verbose == TRUE) message("Computing gene profiles.")
    geneProfiles <- computeGeneProfiles(expressionMatrix = Xoriginal, C = C)
  }

  # End timer and display elapsed time
  end <- Sys.time()
  minutes <- round(difftime(end, start, units = "mins"), digits = 2)
  if (verbose == TRUE) message(paste0("Multivariate curve resolution computed in ", minutes, " minutes."))

  # Create the output list
  outputList <- list("Popt" = P, "Copt" = C, "estimatedMatrix" = Xcalc, "geneProfiles" = geneProfiles)

  return(outputList)
}
