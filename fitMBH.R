fitMBH <- function(x, vars, groups = NULL, nc = 3, ni = 100000, nb = 20000, nt = 20){
  
  library(rjags)
  library(mvtnorm)
  library(mAr)
  library(matrixcalc)
  library(sn)
  
  
  if(any(!vars %in% colnames(x))){stop("Variable names not found in data")}
  
  #get dimensions and covariance matrix
  ndims <- length(vars)
  
  covm <- cov(x[,vars])
  
  if(any(!complete.cases(x))){stop("Missing data found, check all cases are complete before fitting model")}
  
  
  
  
  ##MCMC parameters
  nc = nc # number of chains
  ni = ni # number of iterations
  nb = nb #burnin length 
  nt = nt # thinning parameter
  
  
  #set initial values
  inits <- function() {list(tau = structure(.Data = diag(0.1, ncol = ndims, nrow = ndims), .Dim = c(ndims,ndims)))}
  
  
  
  ########no group model###########
  
  
  if(is.null(groups)){
    
    message("No grouping variable supplied, empirical hypervolume calculated")
    
    nobv <- nrow(x)
    Y2 <- as.matrix(x[,vars])
    cov2 <- covm
    
    #non-nested model
    
    
    sink("mod_empirical.txt")
    cat("model
        {
        #start loop over observations
        for (i in 1:N){
        #each observation has J variables and comes from a multivariate normal described by mu and tau
        Y2[i,1:J] ~ dmnorm(mu[i,1:ndims], tau[ 1:ndims,1:ndims ])
        
        #start loop over variables
        for (j in 1:J){
        #the mean of each variable comes from an independent normal distribution (this allows means to vary between variables)
        mu[i,j] ~ dnorm(0,0.0001)
        }}
        
        
        #the covariance matrix (used to calculate the hypervolume) gets an informative prior based on the covariance matrix calculated from the data
        tau[1 : ndims,1 : ndims] ~ dwish(R[ , ], ndims + 1)
        for (i in 1:ndims){
        for (j in 1:ndims){
        R[i,j] <- cov2[i,j]*(ndims + 1)
        }
        }
        
        }
        ")
    sink()
    
    #specify data
    dat <- list(N = nobv, J = ndims, ndims = ndims,
                Y2 = structure(
                  .Data = as.numeric(Y2),
                  .Dim = c(nobv,ndims)),
                cov2 = structure(
                  .Data = as.numeric(cov2),
                  .Dim = c(ndims,ndims))
    )
    
    #specify parameters to return
    params <- c("tau", "mu")
    
    jm <- jags.model("mod_empirical.txt", data = dat, inits = inits, n.chains = nc, quiet = TRUE)
    
    #burnin
    update(jm, n.iter = nb, progress.bar = "none")
    
    #iterations to use
    
    zm <- coda.samples(jm, variable.names = params, n.iter = ni, thin = nt, progress.bar = "none")
    zj <- jags.samples(jm, variable.names = params, n.iter = ni, thin = nt, progress.bar = "none")
    
    tau <- solve(summary(zj$tau, FUN = mean)$stat)
    mu <- summary(zj$mu, FUN = mean)$stat
    
    ##calculate volume of hypervolume
    
    eig <- eigen(tau)
    
    sf <- qchisq(0.95,df = ndims)
    
    ax <- vector()
    for (k in 1:ndims){
      ax[k] <- sqrt(sf*eig$values[k])
    }
    
    volume <- 2/ndims * (pi^(ndims/2))/factorial((ndims/2)-1) * prod(ax)
    
    
    outlist <- list("means" = mu, "covariance" = tau, "volume" = volume, "Y" = Y2, "dimensions" = vars)
    
  }
  
  
  
 ######## group model #######
  
  
  
  if(!is.null(groups)){
    
    if(!groups %in% colnames(x)){stop("Grouping variable not found in data")}
    
    #nested model
    ngroups <- length(unique(x[,groups]))
    
    #change groups to numeric - ordered alphabetically
    x[,groups] <- as.numeric(as.factor((x[,groups])))
    
    
    nobv <- vector()
      
    for(i in 1:ngroups){
      nobv[i] <- length(x[,groups][x[,groups] == i])
    }
    #observations per group
    
    maxobv <- max(nobv)
      
    #make x into array
      
    xarray <- array(NA, dim = c(maxobv, ngroups, ndims))  
    
    for (j in 1:ndims){
      for (k in 1:ngroups){
      xarray[,k,j] <- c(x[x[,groups]==k, vars[j]], rep(NA, length(xarray[,k,j]) - length(x[x[,groups]==k, vars[j]]))) 
      }
    }
  
    Y <- xarray
    
    
    ##model 1 - with random effect
    
    sink("mod_modelbased.txt")  
    cat("model
        {
        #add new loop over sites 1:K
        for (k in 1:K){
        #loop over observations
        for (i in 1:N){
        
        Y[i,k,1:J] ~ dmnorm(mu[i,k,1:ndims], tau[ 1:ndims,1:ndims ])
        
        #Now allow means to vary by variable and group
        for (j in 1:J){
        mu[i,k,j] ~ dnorm(mu1[k,j], tau.1[j])
        }
        
        } 
        #loop over variables
        for (j in 1:J){
        #each mean comes from an uninformative normal distribution
        mu1[k,j] ~ dnorm(0, 0.0001)
        }
        
        }
        
        #informative prior on covariance matrix for observations
        tau[1 : ndims,1 : ndims] ~ dwish(R[ , ], ndims + 1)
        for (i in 1:ndims){
        for (j in 1:ndims){
        R[i,j] <- covm[i,j]*(ndims + 1)
        }
        }
        
        #uninformative prior on variances for random group effect
        for (j in 1:ndims){    
        tau.1[j] <- 1/pow(sigma.1[j],2)
        sigma.1[j] ~ dunif(0,1)
        }  
        
        }
        ")
    sink()
    
    #specify data required, Y is now a 3d array
    dat <- list(N = maxobv, J = ndims, K = ngroups, ndims = ndims,
                Y = structure(
                  .Data = as.numeric(Y),
                  .Dim = c(maxobv, ngroups ,ndims)),
                covm = structure(
                  .Data = as.numeric(covm),
                  .Dim = c(ndims,ndims))
    )
    
    #specify parameters to return
    params <- c("mu","tau", "tau.1", "mu1")
    
    jm <- jags.model("mod_modelbased.txt", data = dat, inits = inits, n.chains = nc, quiet = TRUE)
    
    #burnin
    update(jm, n.iter = nb, progress.bar = "none")
    
    #iterations to use
    load.module("dic")
    
    zm <- coda.samples(jm, variable.names = params, n.iter = ni, thin = nt, progress.bar = "none")
    zj <- jags.samples(jm, variable.names = params, n.iter = ni, thin = nt, progress.bar = "none")
    
    tau <- solve(summary(zj$tau, FUN = mean)$stat)
    mu <- summary(zj$mu, FUN = mean)$stat
    tau1 <-  1/summary(zj$tau.1, FUN = mean)$stat
    mu1 <- summary(zj$mu1, FUN = mean)$stat
    
    
    ##calculate volume of hypervolume
    
    eig <- eigen(tau)
    
    sf <- qchisq(0.95,df = ndims)
    
    ax <- vector()
    for (k in 1:ndims){
      ax[k] <- sqrt(sf*eig$values[k])
    }
    
    volume <- 2/ndims * (pi^(ndims/2))/factorial((ndims/2)-1) * prod(ax)
    
    
    outlist <- list("means" = mu, "covariance" = tau, "volume" = volume, "group_means" = mu1, "group_variances" = tau1, "Y" = Y, "dimensions" = vars)
    
    
    
    }

  return(outlist)
}
  
  
  
  
  
  

  
  
  
  
  
  
  
  
   
  
  
  
  
  
   
  

  
  
  
  
  
  
  
  
  
  
  

  
  
