plotMBH <- function(x, dims = c(1,2), groupellipses = FALSE, xlim = c(-10,10), ylim = c(-15,15), cols = c("blue", "green", "orange", "purple")){
  
  library(car)
  
  #extract data, means, covariance and variable names from the fitMBH object
  Y <- x$Y
  mu <- x$means
  tau <- x$covariance
  varnames <- x$dimensions

  #define dimensions to display
  d1 <- dims[1]
  d2 <- dims[2]
  
  
  #plot fitted model if no groups are present  
  if(is.null(x$group_means)){
    if(groupellipses == TRUE) {stop("No groups in fitMBH object to show")}
    
    nobs <- nrow(Y)
    
    plot(Y[,d1], Y[,d2], type = "p", xlim = xlim, ylim = ylim, cex.axis = 1.5, cex.lab = 1.5, ylab = varnames[d2], xlab = varnames[d1], pch = 20, col = cols, cex = 2)
    
    ellipse(c(mean(Y[,d1]), mean(Y[,d2])),shape = tau[c(d1,d2), c(d1,d2)],radius = sqrt(2 * qf(.95, 2, nobs)), col = "black", lty = 2)#all data
  
  }
    
    
    
    #plot fitted model if groups are present
    if(!is.null(x$group_means)){
      if(groupellipses == FALSE){
     
        Y2 <- apply(Y, 3, rbind)
        ngroups <- dim(Y)[2]
        nobs <- nrow(Y2)
        
        plot(Y2[,d1], Y2[,d2], type = "n", xlim = xlim, ylim = ylim, cex.axis = 1.5, cex.lab = 1.5, ylab = varnames[d2], xlab = varnames[d1])
        
          for (j in 1:ngroups){
          points(Y[,j,d1], Y[,j,d2], col = cols[j], pch = 20, cex = 2)
          } 
        
        ellipse(c(mean(Y2[,d1]), mean(Y2[,d2])),shape = tau[c(d1,d2), c(d1,d2)],radius = sqrt(2 * qf(.95, 2, nobs)), col = "black", lty = 2)#all data
        
        }
      
      #if groupellipses ==TRUE, plot separate ellipses for each group. This requires calculating a covariance matrix for each group first (note these are not modelled but based on data only)
      if(groupellipses == TRUE){
        
        Y2 <- apply(Y, 3, rbind)
        ngroups <- dim(Y)[2]
        nobs <- dim(Y)[1]
        
        plot(Y2[,d1], Y2[,d2], type = "n", xlim = xlim, ylim = ylim, cex.axis = 1.5, cex.lab = 1.5, ylab = varnames[d2], xlab = varnames[d1])
        
        for (j in 1:ngroups){
          points(Y[,j,d1], Y[,j,d2], col = cols[j], pch = 20, cex = 2)
        }
        
        cov <- list()
        for (j in 1:ngroups){
          cov[[j]] <- cov(Y[,j,])
        }
        
        for (j in 1:ngroups){
           ellipse(c(mean(Y[,j,d1]),mean(Y[,j,d2])),shape = cov[[j]][c(d1,d2), c(d1,d2)],radius = sqrt(2 * qf(.95, 2, nobs)), col = cols[j])
        }
        
      }
      
    }
  
}
    
    
    
    