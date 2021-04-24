#' Function that computes the gene profiles.
#'
#' @author Tobias K. Karakach, Federico Taverna
#'
#' @import data.table
#'
#' @param expressionMatrix The expression matrix organized as (*n* X *m*), *n* is the number of genes and *m* is the number of samples.
#' You should use the expressionMatrix output of prepareData().
#'
#' @param C The C matrix. You should use the Copt output of multivariateCurveResolution().
#'
#' @return The gene profiles data.frame.
computeGeneProfiles <- function(expressionMatrix, C)
{
  # Define the cosine similarity function
  cos.sim <- function(x, y)
  {
    A <- x
    B <- y
    return(sum(A * B) / sqrt(sum( A ^ 2) * sum( B ^ 2)))
  }

  # Assign the arguments to new variable names
  Xorg <- t(expressionMatrix)
  ncomp <- ncol(C)
  rm(expressionMatrix)

  # Calculate the cosine similarity
  ngenes <- ncol(Xorg)
  angles <- matrix(0, ncomp, ngenes)
  colnames(angles) <- colnames(Xorg)

  for (j in 1: ncomp)
  {
    for (i in 1:ngenes)
    {
      angles[j, i] <- cos.sim(C[, j], Xorg[, i])
    }
  }

  # Omit rows with NAs
  angles <- angles[, colSums(is.na(angles)) == 0]

  # Generate the gene profiles
  geneProfiles <- data.frame((matrix(vector(), dim(angles)[2], ncomp, dimnames = list(c(), colnames(C)))))

  for (i in 1:ncomp)
  {
    k <- sort(angles[i, ], decreasing = TRUE, na.last = TRUE)
    geneProfiles[, i] <- names(k)
  }

  return (geneProfiles)
}
