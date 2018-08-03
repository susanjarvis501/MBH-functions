##run MBH functions

source("simulate_dataMBH.R")
source("fitMBH.R")
source("plotMBH.R")
source("overlapMBH.R")
source("inclMBH.R")


#simulate data 

data1 <- simulate_dataMBH(nobs = 20, ndims =3, sdgrp = 2, ngroups = 4, variances = "fixed",  returntruecov = TRUE)

#fit empirical hypervolume ignoring group structure

hvE <- fitMBH(data1$data, vars = c("V1", "V2", "V3"), groups = NULL)

#fit model-based hypervolume

hvMB <- fitMBH(data1$data, vars = c("V1", "V2", "V3"), groups = "Group")

#plot data with group ellipses

plotMBH(hvMB, dims = c(1,2), groupellipses = TRUE, xlim = c(-10, 12))

#plot empirical hypervolume

plotMBH(hvE, dims = c(1,2))

#plot model-based hypervolume

plotMBH(hvMB, dims = c(1,2))

#add true covariance matrix for comparison

ellipse(c(mean(data1$data[,1]), mean(data1$data[,2])),shape = data1$truecov[c(1,2), c(1,2)],radius = sqrt(2 * qf(.95, 2, nrow(data1$data))), col = "red", lty = 2)


#simulate a second dataset and calculate overlap

data2 <- simulate_dataMBH(means = c(3,3,3), nobs = 10, ndims =3, sdgrp = 3, ngroups = 5, variances = "fixed",  returntruecov = TRUE)

#fit model-based hypervolume

hvMB2 <- fitMBH(data2$data, vars = c("V1", "V2", "V3"), groups = "Group")

#calculate overlap and plot (caution - this is quite slow)

overlapMBH(hvMB, hvMB2, plot = TRUE)

#create some new data and calculate inclusion probabilities

new_data <- simulate_dataMBH(nobs = 10, ngroups = 1, ndims = 3, means = c(5,5,5))
new_data

#plot new data (note use overlap = FALSE argument to only plot overlap and not calculate statistic)
overlapMBH(hvMB, hvMB2, overlap = FALSE, plot = TRUE, dims = c(1,2))
points(new_data$data[,1], new_data$data[,2], pch = 20)

#test inclusion in hypervolume 1
inclMBH(hvMB, new_data$data)

#test inclusion in hypervolume 2
inclMBH(hvMB2, new_data$data)




