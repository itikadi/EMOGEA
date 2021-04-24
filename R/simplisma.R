#' Function to obtain the initial values for P or C.
#' Ref: W. Windig, J. Guilment, Anal. Chem., 63, 1425 (1991)
#'
#' @author Tobias K. Karakach, Federico Taverna
#'
#' @import data.table
#' @importFrom stats sd
#'
#' @param expressionMatrix The expression matrix organized as (*n* X *m*), *n* is the number of genes and *m* is the number of samples.
#' You should use the expressionMatrix output of prepareData().
#'
#' @param numberOfComponents The number of number of components (profiles) to be extracted.
#'
#' @param noiseFactor The noise factor.
#'
#' @return A list which contains the initial estimates for P (Pinit) and C (Cinit).
simplisma <- function(expressionMatrix, numberOfComponents, noiseFactor = 0.009)
{
  # Assign the arguments to new variable names
  data <- t(expressionMatrix)
  comp <- numberOfComponents
  rm(expressionMatrix, numberOfComponents)

  # Algorithm initialization
  data[is.na(data)] <- 0
  r <- dim(data)[1]
  c <- dim(data)[2]
  alpha <- noiseFactor          # Noise factor
  ipure <- matrix(1, 1, comp)   # Lists purest variables
  var.sim <- matrix(0, 1, c)    # Lists similarity for each variable
  coo <- matrix(0, comp, comp)  # Correlation around the origin matrix
  lamda_i <- sqrt((apply(data, 2, sum) ^ 2) / length(data[, 1])) # Length of variables (Eqn 1)

  # Pure variable selection
  mu_i <- apply(data, 2, mean)    # Mean value for each channel (Eqn 3)
  mean.data <- apply(data,2,mean) # Mean value for each channel (Eqn 3)
  sigma_i <- apply(data, 2, sd)   # Standard deviation of variable (Eqn 4)
  max.mean.alpha <- max(mu_i * alpha)         # Fudge factor for baseline noise
  purity <- sigma_i / (mu_i + max.mean.alpha) # Eqn 6
  d_norm <- sqrt((mu_i ^ 2 + (sigma_i + max.mean.alpha) ^ 2))
  D <- sweep(data, 2, d_norm, "/")            # Eqn 8
  s <- matrix(0, comp, c)

  # Start loop
  for (i in 1:comp)
  {
    pmax <- -1e99   # Stores value of maximum purity channel

    # Calculate length of vector from each channel to an origin
    # 1st variable: choose longest vector from center of plane
    # Others: choose longest vector away from previous chosen pure variables
    for (j in 1:c)
    {
      for (k in 1:i)
      {
        z <- sum(data[, j] * data[, ipure[k]]) / (d_norm[j] * d_norm[ipure[k]])
        coo[k, i] <- z
        coo[i, k] <- z
      }

      coo[i,i] <- sum(data[, j] ^ 2) / (d_norm[j] ^ 2)

      # Determinant of coo gives similarity between columns; lower det is similar
      if (coo[i,i] > 1e-6)
      {
        var.sim[j] <- det(as.matrix(coo[1:i, 1:i]))
      }
      else
      {
        var.sim[j] <- 0 # Error trap
      }

      # Keep track of purest channel so far (highest product)
      if (var.sim[j] * purity[j] > pmax)
      {
        ipure[i] <- j		# If channel j is purest, save it
        pmax <- var.sim[j] * purity[j]
      }

      s[i,j] <- sigma_i[j] * var.sim[j]
    }

    w1 <- lamda_i ^ 2 / (mean.data ^ 2 + (sigma_i + alpha) ^ 2)
    s1 <- sigma_i * w1

    # Once the pure variable is chosen, coo is set for next component
    for (k in 1:i)
    {
      z <- sum(data[, ipure[i]] * data[, ipure[k]]) / (d_norm[ipure[i]] * d_norm[ipure[k]])
      coo[k,i] <- z
      coo[i,k] <- z
    }

    coo[i,i] <- sum(data[, ipure[i]] ^ 2) / (d_norm[ipure[i]] ^ 2)
  }

  # Rs <- 100 * (colSums(t(s) / sum(s1)))
  # Rr <- matrix(0,1,comp)
  # for (i in 1:comp-1) Rr[i] <- Rs[i] / Rs[i+1]

  # Compute C and P
  purenames <- colnames(data)[ipure]
  C	<- data[, ipure]
  colnames(C) <- purenames
  rownames(C) <- rownames(data)
  S	<- t(data) %*% C %*% solve(t(C) %*% C)
  P <- t(S)

  # Create the output list
  outputList <- list("Cinit" = C, "Pinit" = P)

  return(outputList)
}




