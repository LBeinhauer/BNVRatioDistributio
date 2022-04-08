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
