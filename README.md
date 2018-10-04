# MBH functions

Functions to create model-based hypervolumes.

## Examples of how to run each function can be found in [run MBH functions.R](https://github.com/susanjarvis501/MBH-functions/blob/master/run%20MBH%20functions.R)



**Required packages**

rjags  
mvtnorm  
mAr  
matrixcalc  
sn  
car  
    


## [simulate_dataMBH](https://github.com/susanjarvis501/MBH-functions/blob/master/simulate_dataMBH.R)

`simulate_dataMBH(nobs = 10, ndims = 3, ngroups = 4, sdgrp = 2, sdobs = NULL, variances = "fixed", vardiff = 0.5, means = rep(0,ndims), returntruecov = FALSE)`

**Arguments**

`nobs` - Number of observations to simulate per group  
`ndims` - Number of variables to simulate  
`ngroups` - Number of groups to simulate  
`sdgrp` - Standard deviation of group means  
`sdobs` - Standard deviation of observations within groups. Either a single value or vector with length equal to ndims. If NULL (default) a vector of length equal to ndims is generated from a Uniform distribution with bounds 1 and 10. Note that if sdobs = NULL the vector generated is on the variance scale, if sdobs is provided then it should be a standard deviation
`variances` - Either "fixed" to constrain within-group variances to be the same (the default), or "variable" so within-group variances can vary  
`vardiff` - Standard deviation of within-group variance differences  
`means` - Means of each variable to simulate  
`returntruecov` - Logical. Return true covariance matrix as well as data?  

**Value**

`data` - simulated data  
`truecov` - covariance matrix used to simulate data  

**Details**

Returns a list


## [fitMBH](https://github.com/susanjarvis501/MBH-functions/blob/master/fitMBH.R)

`fitMBH(x, vars, groups = NULL, nc = 3, ni = 100000, nb = 20000, nt = 20)`

**Arguments**

`x` - Data to fit hypervolume to  
`vars` - Names of variables in x to use for hypervolume construction  
`groups` - Name of the grouping variable in x. Use NULL if no groups present or to  ignore grouping structure and fit an empirical hypervolume  
`nc` - Number of MCMC chains  
`ni` - Number of MCMC iterations (default 100000)  
`nb` - Length of burnin    
`nt` - Thinning paramter    

**Value**

`means` - Estimated means of each variable  
`covariance` - Estimated covariance structure  
`volume` - Estimated hypervolume size  
`group_means` - Estimated group means  
`group_variances` - Estimated between-group variances for each variable  


## [plotMBH](https://github.com/susanjarvis501/MBH-functions/blob/master/plotMBH.R)

`plotMBH(x, dims = c(1,2), groupellipses = FALSE, xlim = c(-10,10), ylim = c(-15,15), cols = c("blue", "green", "orange", "purple"))`

**Arguments**

`x` - Fitted MBH model  
`dims` - Which two dimensions to plot  
`groupellipses` - Logical. Plot ellipses for each group?  
`xlim` - xlim of plot  
`ylim` - ylim of plot  
`cols` - Colours to plot groups  


## [overlapMBH](https://github.com/susanjarvis501/MBH-functions/blob/master/overlapMBH.R)

`overlapMBH(hv1, hv2, overlap = TRUE, plot = TRUE, dims = c(1,2), col1 = "black", col2 = "blue")`

**Arguments**

`hv1` - Fitted MBH model  
`hv2` - Fitted MBH model   
`overlap` - Logical. Do you want to calculate overlap? This can be very slow  
`plot` - Logical. Do you want to plot overlap?  
`dims` - Dimensions to plot  
`col1` - Colour to use for first hypervolume  
`col2` - Colour to use for second hypervolume  

**Details**

Utilises a simulation based approach to calculate overlap by simulating a number of points from each hypervolume. Returns an overlap statistic defined by total number of points shared divided by total number of points simulated. The density of points in each hypervolume is kept constant. Can be very slow for large hypervolumes.



## [inclMBH](https://github.com/susanjarvis501/MBH-functions/blob/master/inclMBH.R)

`inclMBH(hv, newdat)`

**Arguments**

`hv` - Fitted MBH model  
`newdat` - New data to test probability of inclusion in the hypervolume

**Details**

Returns newdat with an additional column corresponding to probability of inclusion in the hypervolume





























