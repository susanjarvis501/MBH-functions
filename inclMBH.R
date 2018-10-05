inclMBH <- function(hv, newdat){
  
  library(car)
  library(mvtnorm)
  
  #extract volume
  
  vol1 <- hv$volume
  
  #simulate points from each hypervolume
  
  #calculate means
  
  if(is.null(hv$group_means)){
    mean1 <- colMeans(hv$means)
    
  }else{mean1 <- colMeans(apply(hv$means, 3, rbind))}
  
  
  #extract covariance
  
  cov1 <- hv$covariance
  
  #extract variable names
  varnames1 <- hv$dimensions
  
  #extract correct vars from newdat
  newvars <- newdat[,varnames1]
  
  #generate random points from hypervolume
  pnts_hv <- rmvnorm(round(vol1), mean1, cov1, method = "eigen")
  
  
  #test inclusion of new points
  totestall <- newvars
  prob <- vector()
  mean.test.p <- vector()
  for(k in 1:nrow(totestall)){
    totest <- as.numeric(totestall[k,])
    test.p <- vector()
    tau <- cov1
    #colnames(tau) <- rownames(tau)
    mu <- mean1
    #test new point against distribution
    prob <- min(pmvnorm(upper = totest,sigma = tau, mean = mu),pmvnorm(lower = totest,sigma = tau, mean = mu)*2)
    
    #simulate 99 draws from multivariate dist (for speed)
    rsims <- pnts_hv[sample(nrow(pnts_hv), 99, replace = FALSE), ]
    #calculate p values for each simulation
    sim.prob <- vector()
    for (j in 1:nrow(rsims)){
      sim.prob[j] <- min(pmvnorm(upper = rsims[j,],sigma = tau, mean = mu),pmvnorm(lower = rsims[j,],sigma = tau, mean = mu*2))
    }
    #calc probability of inclusion
    all.prob <- c(sim.prob,prob)
    prob.df <- ecdf(all.prob)
    #plot
    #plot(prob.df); abline(v=prob[i])
    test.p <- prob.df(prob)
    mean.test.p[k] <- mean(test.p)
  }
  
  p.out <- cbind(totestall, mean.test.p)
  
  return(p.out)
  
  
  
}