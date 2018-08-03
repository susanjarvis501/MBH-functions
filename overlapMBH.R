overlapMBH <- function(hv1, hv2, overlap = TRUE, plot = TRUE, dims = c(1,2), col1 = "black", col2 = "blue"){
  
  library(car)
  library(mvtnorm)
  
  #extract volumes
  
  vol1 <- hv1$volume
  vol2 <- hv2$volume
  
  #simulate points from each hypervolume
  
  #calculate means
  
  if(is.null(hv1$group_means)){
    mean1 <- colMeans(hv1$means)
    
  }else{mean1 <- colMeans(apply(hv1$means, 3, rbind))}
  
  if(is.null(hv2$group_means)){
    mean2 <- colMeans(hv2$means)
    
  }else{mean2 <- colMeans(apply(hv2$means, 3, rbind))}
  
  
  #extract covariances
  
  cov1 <- hv1$covariance
  
  cov2 <- hv2$covariance
  
  #extract variable names
  varnames1 <- hv1$dimensions
  varnames2 <- hv2$dimensions
  
  if(!identical(varnames1, varnames2)) {stop ("Different variables in each hypervolume or variables in different orders")}
  
  
  if(overlap == TRUE){

  #simulate random points from each distribution
  
  pnts_hv1 <- rmvnorm(round(vol1), mean1, cov1, method = "eigen")
  pnts_hv2 <- rmvnorm(round(vol2), mean2, cov2, method = "eigen")
  
  #check if points are in each distribution
  
  #hv1 points in hv2 
  
  totestall <- pnts_hv1
  prob <- vector()
  mean.test.p <- vector()
  for(k in 1:nrow(totestall)){
    totest <- as.numeric(totestall[k,])
    test.p <- vector()
    tau <- cov2
    #colnames(tau) <- rownames(tau)
    mu <- mean2
    #test new point against distribution
    prob <- min(pmvnorm(upper = totest,sigma = tau, mean = mu),pmvnorm(lower = totest,sigma = tau, mean = mu)*2)
    
    #simulate 999 draws from multivariate dist
    rsims <- pnts_hv2[sample(nrow(pnts_hv2), 99, replace = FALSE), ]
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
  
  no_hv1_in_hv2 <- nrow(p.out[p.out[,4]>0.05,])
  
  
  
  #hv2 points in hv1 
  
  totestall <- pnts_hv2
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
    
    #simulate 999 draws from multivariate dist
    rsims <- pnts_hv1[sample(nrow(pnts_hv1), 99, replace = FALSE), ]
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
  
  no_hv2_in_hv1 <- nrow(p.out[p.out[,4]>0.05,])
  
  overlap_12 <- (no_hv1_in_hv2 + no_hv2_in_hv1)/(vol1 + vol2)
  
  return(overlap_12)
  
  }
  
  if(plot == TRUE){
    
    #define dimensions to display
    d1 <- dims[1]
    d2 <- dims[2]
    
    #format data
    if(length(dim(hv1$Y)) > 2){
      hv1Y <- apply(hv1$Y, 3, rbind)
    } else {hv1Y <- hv1$Y}
    
    if(length(dim(hv2$Y)) > 2){
      hv2Y <- apply(hv2$Y, 3, rbind)
    } else {hv2Y <- hv2$Y}
    
    nobs1 <- nrow(hv1Y)
    nobs2 <- nrow(hv2Y)
    
    minX <- min(min(hv1Y[,d1]), min(hv2Y[,d1])) - 0.2*(range(hv1Y[,d1], hv2Y[,d1])[2]-range(hv1Y[,d1], hv2Y[,d1])[1])
    minY <- min(min(hv1Y[,d2]), min(hv2Y[,d2])) - 0.2*(range(hv1Y[,d2], hv2Y[,d2])[2]-range(hv1Y[,d2], hv2Y[,d2])[1])
    maxX <- max(max(hv1Y[,d1]), max(hv2Y[,d1])) + 0.2*(range(hv1Y[,d1], hv2Y[,d1])[2]-range(hv1Y[,d1], hv2Y[,d1])[1])
    maxY <- max(max(hv1Y[,d2]), max(hv2Y[,d2])) + 0.2*(range(hv1Y[,d2], hv2Y[,d2])[2]-range(hv1Y[,d2], hv2Y[,d2])[1])
    
      
    plot(hv1Y[,d1], hv1Y[,d2], type = "n", xlim = c(minX, maxX), ylim = c(minY, maxY), cex.axis = 1.5, cex.lab = 1.5, ylab = varnames1[d2], xlab = varnames1[d1])
      

   ellipse(c(mean1[1], mean1[2]),shape = cov1[c(d1,d2), c(d1,d2)],radius = sqrt(2 * qf(.95, 2, nobs1)), col = col1, lty = 2, fill = TRUE)
   
   ellipse(c(mean2[1], mean2[2]),shape = cov2[c(d1,d2), c(d1,d2)],radius = sqrt(2 * qf(.95, 2, nobs2)), col = col2, lty = 2, fill = TRUE)
   
   
    }
  
  
    
  }
  
  
