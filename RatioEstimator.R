# Estimating Mean and Variance of Ratio of bivariate gaussian variables #

## van Kempen LJ van Vliet (1999)




# unbiased estimator r3 of R = x/Y

r_3 <- function(x_bar, y_bar, n, mu_x, mu_y, var_y, cov_xy){
  (x_bar / y_bar) - (1/n) * (((mu_x/mu_y^3) * var_y) - ((cov_xy)/(mu_y^2)))
}


# variance of unbiased estimator r3 (?)

var_r3 <- function(n, mu_x, mu_y, var_x, var_y, cov_xy){
  (1/n) * ((var_x/mu_y^2) + ((mu_x^2 * var_y)/mu_y^4) - ((2 * mu_x * cov_xy)/mu_y^3))
}






# Testing estimator using simulated data

library(MASS)


# assemble correlation matrix
cor <- matrix(c(1, .5, .5, 1), nrow = 2, byrow = T)

# assemble covariance matrix
covmat <- cor %*% diag(c(3^2, 3^2))

# generate data, using covariance matrix and means of 10 for both variables
dat <- mvrnorm(n = 1000, mu = c(10, 10), Sigma = covmat)

# "true" values of mean of ratio
mean(dat[,1])/mean(dat[,2]) # mean_x / mean_y
mean(dat[,1]/dat[,2])       # mean_r
var(dat[,1]/dat[,2])        # var_r

# estimated ratio
r_3(mean(dat[,1]), mean(dat[,2]), 1000, 10, 10, var(dat[,2]), cov(dat[,1], dat[,2]))

# variance of ratio estimator
var_r3(1000, 10, 10, var(dat[,1]), var(dat[,2]), cov(dat[,1], dat[,2]))




var_ratio <- function(var_x, var_y, mu_x, mu_y, cov_xy){
  var_x/mu_y^2 + (mu_x^2 * var_y)/mu_y^4 - (2 * mu_x * cov_xy)/mu_y^3
}


sim_ratio <- function(var_x, var_y, mu_x, mu_y, cov_xy, n = 1000){
  cm <- matrix(cov_xy, nrow = 2, ncol = 2)
  diag(cm) <- c(var_x, var_y)
  
  data <- MASS::mvrnorm(n, mu = c(mu_x, mu_y), Sigma = cm)
  
  obsm <- mean(data[,1]/data[,2])
  estm <- r_3(mean(data[,1]), mean(data[,2]), n, mu_x, mu_y, var_y, cov_xy)
  
  obsv <- var((data[,1]/data[,2]))
  estv <- var_ratio(var_x, var_y, mu_x, mu_y, cov_xy)
  
  df <- data.frame(var_x = var_x,
                   var_y = var_y,
                   mu_x = mu_x,
                   mu_y = mu_y,
                   cov_xy = cov_xy,
                   obsm = obsm,
                   estm = estm,
                   obsv = obsv,
                   estv = estv)
  
  df
}


simmat <- matrix(NA, nrow = 1000, ncol = 10)

vx <- seq(.01, 10, by = .01)

for(i in 1:1000){
  sr <- sim_ratio(var_x = vx[i], var_y = 3, mu_x = 10, mu_y = 10, cov_xy = 0, n = 10000)
  simmat[i,1] <- sr$obsv
  simmat[i,2] <- sr$estv
}

vy <- seq(.01, 10, length.out = 1000)

for(i in 1:1000){
  sr <- sim_ratio(var_x = 3, var_y = vy[i], mu_x = 10, mu_y = 10, cov_xy = 0, n = 10000)
  simmat[i,3] <- sr$obsv
  simmat[i,4] <- sr$estv
}


mux <- seq(-10, 10, length.out = 1000)

for(i in 1:1000){
  sr <- sim_ratio(var_x = 3, var_y = 3, mu_x = mux[i], mu_y = 10, cov_xy = 0, n = 10000)
  simmat[i,5] <- sr$obsv
  simmat[i,6] <- sr$estv
}


muy <- seq(.01, 10, length.out = 1000) # variance not negative!

for(i in 1:1000){
  sr <- sim_ratio(var_x = 3, var_y = 3, mu_x = 10, mu_y = muy[i], cov_xy = 0, n = 10000)
  simmat[i,7] <- sr$obsv
  simmat[i,8] <- sr$estv
}


covxy <- seq(-3, 3, length.out = 1000)

for(i in 1:1000){
  sr <- sim_ratio(var_x = 3, var_y = 3, mu_x = 10, mu_y = 10, cov_xy = covxy[i], n = 10000)
  simmat[i,9] <- sr$obsv
  simmat[i,10] <- sr$estv
}




plot(simmat[,1], simmat[,2], xlim = c(.03, .16)) # estimate about .8 of size of observed 
plot(simmat[,3], simmat[,4], xlim = c(.0255, .28)) # non-linear, if variance in y is multiple of var in x
plot(simmat[,5], simmat[,6]) # estimate about .9 of size of observed 
plot(simmat[,7], simmat[,8], xlim = c(0, 1), ylim = c(0, .3))
plot(simmat[,9], simmat[,10]) # estimate about .8 of size of observed 

