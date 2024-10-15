#===============================================================================
# Author: [Rana Bamdadi]
# Description: This script estimates the parameters of a Geometric time series 
#              model using the Whittle likelihood function and evaluates model 
#              performance with various error metrics.
#===============================================================================
############################################
################## G-MINBL(1,0,1,1) model 
################## Geometric Modified Integer-Valued Bilinear
############################################
# Clear the workspace
rm(list = ls(all = TRUE))

# Load required libraries
library(TSA)
library(stats)
library(rootSolve)
library(readxl)

#===============================================================================
# Data Loading
#===============================================================================
# Load data from the Excel file
# Ensure that the file path is correct for your local system
my_data <- read_excel("/Users/Documents/Data/RCI_offencebymonth.xlsm")
x <- as.numeric(my_data[3201, 4:354])  # Extract the relevant data

N <- length(x)  # Length of the time series

#===============================================================================
# Autocovariance Calculation
#===============================================================================
# Calculate autocovariances at lags 0, 1, and 2
g0 <- acf(x, lag.max = 50, type = "covariance", plot = FALSE)$acf[1]
g1 <- acf(x, lag.max = 50, type = "covariance", plot = FALSE)$acf[2]
g2 <- acf(x, lag.max = 50, type = "covariance", plot = FALSE)$acf[3]

#===============================================================================
# Initial Parameter Estimates
#===============================================================================
fn <- function(a) 
  c(
    F1=a[1]*a[2]*((1-a[3])/(a[3]^2))+((1-a[3])/a[3])-mean(x)*(1-(a[1]+a[1]*a[2]*((1-a[3])/a[3]))),
    F2=(a[1]+a[1]*a[2]*((1-a[3])/a[3]))-g2/g1,
    F3=a[1]*a[2]*(2*(((1-a[3])/a[3])^3)+3*(((1-a[3])/a[3])^2)+((1-a[3])/a[3])+((1-a[3])/(a[3]^2))*mean(x))-(g1-(a[1]+a[1]*a[2]*((1-a[3])/a[3]))*g0)
  )
ss <- multiroot(f = fn, start =c(0.315963,-0.30969,0.634464))

aahatx=ss$root[1]
bbhatx=ss$root[2]
muhatx=ss$root[3]

#===============================================================================
# Periodogram Calculation
#===============================================================================
# Calculate periodogram and associated frequencies
h <- (N - 1) / 2  # Half-length of the series
Iprio <- periodogram(x, log = 'no', plot = FALSE)$spec[1:h]  # Periodogram
w <- (2 * pi * (1:h)) / (N - 1)  # Frequencies

#===============================================================================
# Whittle Likelihood Function
#===============================================================================
# Define the Whittle likelihood function for parameter estimation
whittle_likelihood <- function(l) {
  ahatx <- l[1]
  bhatx <- l[2]
  muhatx <- (1-l[3])/l[3]
  log_likelihood <- 0
  
  for (k in 1:h) {
    # Spectral density function based on Whittle likelihood
    f_X_wk <- (1 / (2 * pi)) *(g0 + 2 * g1 * ((cos(w[k]) - (ahatx + (ahatx * bhatx * muhatx))) /
                                                (1 + (((ahatx + (ahatx * bhatx * muhatx))^2) - 2 * ((ahatx + (ahatx * bhatx * muhatx)) * cos(w[k]))))))
    log_likelihood <- log_likelihood + (log(f_X_wk) + (Iprio[k] / f_X_wk))
  }
  return(log_likelihood / N)  # Return normalized log-likelihood
}

#===============================================================================
# Parameter Estimation using Non-Linear Minimization (NLM)
#===============================================================================
# Initial parameter vector
initial_params <- c(aahatx,bbhatx,muhatx)

# Estimate parameters using nlm
est <- nlm(whittle_likelihood, initial_params)

# Extract estimated parameters
a <- est$estimate[1]
b <- est$estimate[2]
mu_eps <- est$estimate[3]

#===============================================================================
# Residuals and Forecast Computation
#===============================================================================
# Initialize arrays for estimated values (xhat) and residuals (e)
xhat <- rep(NA, N)
e <- rep(NA, N)
e[1] <- 1  # Set the initial residual

# Compute residuals and estimated values iteratively
for (i in 2:N) {
  e[i] <- round(max(0, (x[i] - (a * x[i - 1]) - (a * b * x[i - 1] * e[i - 1]))))
}

xhat[1]=(1-mu_eps)/mu_eps
for(j in 2:N){
  xhat[j]=(1-mu_eps)/mu_eps+(a*x[j-1])+(a*b*x[j-1]*e[j-1])
}

#===============================================================================
# Model Performance Evaluation
#===============================================================================
# Calculate error metrics: Root Mean Square Error, Mean Absolute Error, etc.
RMS <- sqrt(mean(((abs(x - xhat))^2)))
MAE <- mean(abs(x - xhat))
MdAE <- median(abs(x - xhat))
MSE <- mean((x - xhat)^2)

#===============================================================================
# In-Sample Forecasting Error (Last 10 Periods)
#===============================================================================
# Calculate forecast errors for the last 10 periods
FMSE10 <- mean((x[342:351] - xhat[342:351])^2)   # Forecast MSE for 10 periods
FMAE10 <- mean(abs(x[342:351] - xhat[342:351]))  # Forecast MAE for 10 periods

#===============================================================================
# End of Script
#===============================================================================