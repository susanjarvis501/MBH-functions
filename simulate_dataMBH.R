simulate_dataMBH <- function(nobs = 10, ndims = 3, ngroups = 4, sdgrp = 2, variances = "fixed", vardiff = 0.5, means = rep(0,ndims), returntruecov = FALSE){
  library(mvtnorm)
  
  #create a random covariance matrix
  n <- ndims
  p <- qr.Q(qr(matrix(rnorm(n^2), n)))
  truecov <- crossprod(p, p*(runif(ndims,1,10)))
  
  #create an empty array for the data
  Y <- array(NA, dim = c(nobs,ngroups,ndims))
  
  #if within-group variances are fixed, simulate the required number of observations from a multivariate normal with the true covariance matrix and a random mean vector (with standard deviation defined by between-group variance parameter)
  if (variances == "fixed"){
    
    for(j in 1:ngroups){
      Y[,j,] <- rmvnorm(nobs, mean = rnorm(ndims,means,sdgrp), sigma = truecov)
    }
  }
  
  #if within-group variances are fixed, add some variation to the within-group variances, then simulate the required number of observations as above
  if (variances == "variable"){
    
    for(j in 1:ngroups){
      truecor <- cov2cor(truecov)
      sds <- sqrt(diag(truecov))
      sds1 <- sqrt(diag(truecov)) + runif(3,-vardiff, vardiff)
      newcov <- sweep(sweep(truecor, 1L, sds1, "*"), 2L, sds1, "*")
      Y[,j,] <- rmvnorm(nobs, mean = rnorm(ndims,0,sdgrp), sigma = newcov)
    }
  }
  
  #condense the array to a dataframe
  Ydf <- as.data.frame(apply(Y,3L,c))
  
  #add a field to denote group membership
  Ydf$Group <- rep(1:ngroups, each = nobs)
  
  #return the dataframe and (if specified) the true covariance matrix
  if(returntruecov == TRUE){
    outlist <- list("data" = Ydf, "truecov" = truecov)
    return(outlist)
  } else{
    outlist <- list("data" = Ydf)
    return(outlist)
  }
}

#dat1 <- simulate_dataMBH()
#dat1
