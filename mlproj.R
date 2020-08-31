#
#' This is a function for Maximum Likelihood PCA with equal row correlated errors
#' x is the data matrix organized as (mxn), where m is the number of samples and n is the number of features
#' xcov is the inverse of the measurement errors whose dimensions are (n x n), and
#' ncomp is the number of principal components estimated to be required to reconstruct x without loss of generality
#
#' @param x data matrix
#' @param xcov covariance matrix
#' @param ncomp number of components
#' @return A matrix of the infile
#' @export
#'
mlproj<-function(x,xcov,ncomp)
{
  m<-dim(x)[1]
  n<-dim(x)[2]
  pc <- ncomp
  flg <- 1
  maxiter <- 2000
  lamda <- 1e-10
  iter <- 0
  #
  tmp=svd(x)
  u<-tmp$u
  v<-tmp$v
  s<-diag(tmp$d)
  #
  V <- tmp$v[,1:pc]
  while (flg ==1 || iter<maxiter)
  {
    iter <- iter+1
    Xhat <- x%*%((Xcov))%*%V%*%solve(t(V)%*%(Xcov)%*%V)%*%t(V)
    S1  <- sum(diag(((x-Xhat)%*%(Xcov)%*%t((x-Xhat)))))
    #
    tmp <-svd(Xhat)
    Xhat <- tmp$u[,1:pc]%*%t(tmp$u[,1:pc])%*%x
    S2   <- sum(diag(((x-Xhat)%*%(Xcov)%*%t((x-Xhat)))))
    tmp<-svd(Xhat)
    U <- tmp$u[,1:pc]
    S <- diag(tmp$d[1:pc])
    V <- tmp$v[,1:pc]
    #
    sobj <- abs((S1-S2)/S2)
    if(sobj < lamda){
      flg = 0
      break
    }
    print(paste0("Iteration number: ", iter, "; Converging to: ",lamda, "; Currently at = ", sobj))
    output <- list("U" = U,
                   "S" = S,
                   "V"=V)
  }
  return(output)
}
