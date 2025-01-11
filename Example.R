library(Rcpp)
library(Rfast)
library(zoo)
Rcpp::sourceCpp("the algorthims.cpp")


# simulate an example 
set.seed(111)
h <- c( rvonmises(1000 , -pi/2, 10) , rvonmises(1000 , pi/2, 10) )
# Calculate rho , kappa_rho , kappa_mu. 

## Using secussive difference 
rho <- median(1-cos(diff(h,2))) / median(1-cos(diff(h,1))) -1 
kappa_rho <- 0.4592 / ( (1+rho)*(1-sqrt( median(sin(diff(h)))^2 + median(cos(diff(h)))^2 )))
kappa_mu <- (1-rho^2)*kappa_rho

# Using non-parameteric method

# Compute rolling medians
ms <- rollmedian(sin(h), window, na.pad = TRUE, align = "center")
mc <- rollmedian(cos(h), window, na.pad = TRUE, align = "center")

# Compute moving median 
med_head <- atan2(ms, mc)

# residuals of moving median filter
errors <- na.omit((h - med_head) %% (2 * pi))
sin_errors <- sin(errors)
cos_errors <- cos(errors)
n <- length(errors)

# Compute rho (autocorrelation)
rho <- sum(sin_errors[-1] * sin_errors[-n]) / 
        sqrt(sum(sin_errors[-1]^2) * sum(sin_errors[-n]^2))

# Calculate whitened errors
whitened_errors <- (errors[-1] - atan2(rho * sin_errors[-n], 1 - rho + rho * cos_errors[-n])) %% (2 * pi)

# Estimate kappa for whitened_h and errors
kappa_rho <- Rfast::vm.mle(whitened_errors)$param[2]
kappa_mu  <- Rfast::vm.mle(errors)$param[2]




# fit the IID model 
penalty <-  -2*log(length(h))/kappa_mu   * (1+rho)/(1-rho)
tps <- AFPOP(sin(h), cos(h), penalty)
sort(tps)

