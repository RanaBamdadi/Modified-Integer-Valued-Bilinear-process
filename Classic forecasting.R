#===============================================================================
################## Forecasting study
# The code in this file reproduces the Classic forecasting results 
#===============================================================================
############################################
################## PL-MINBL(1,0,1,1) model 
################## Poisson_lindly Modified Integer-Valued Bilinear
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
    F1=(a[1]*a[2]*(((2*((a[3] + 2)/(a[3]*(a[3] + 1)))) + (2*((a[3] + 3)/(((a[3])^2)*(a[3] + 1))))) - (((a[3] + 2)/(a[3]*(a[3] + 1)))^2)))+((a[3] + 2)/(a[3]*(a[3] + 1)))-mean (x)*(1 - (a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))),
    F2=(a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))-g2/g1,
    F3=(a[1]*a[2]*((6*((a[3] +4)/(((a[3])^3)*(a[3]+ 1)))) + (12*((a[3] +3)/(((a[3])^2)*(a[3]+ 1)))) + (5*(a[3]+2)/(a[3]*(a[3]+ 1)))))+
      (a[1]*a[2]*mean (x)*(((2*(a[3]+2))/(a[3]*(a[3]+1)))+((2*(a[3]+3)/((a[3]^2)*(a[3]+1))))-(((a[3]+2)/(a[3]*(a[3]+1)))^2)))+
      (-3*a[1]*a[2]*((a[3] +2)/(a[3]*(a[3] +1)))*(((2*(a[3]+2))/(a[3]*(a[3]+1)))+((2*(a[3]+3)/((a[3]^2)*(a[3]+1))))-(((a[3]+2)/(a[3]*(a[3]+1)))^2)))-
      (a[1]*a[2]*(((a[3] + 2)/(a[3]*(a[3] + 1)))^3)) - (g1 - ((a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))*g0)) 
  )
ss <- multiroot(f = fn, start =c(0.284446,-0.159125,2.3623))

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
  muhatx <- (l[3]+2)/(l[3]*(l[3]+1))
  log_likelihood <- 0
  
  for (k in 1:h) {
    # Spectral density function based on Whittle likelihood
    f_X_wk <- (1 / (2 * pi)) *(g0 + 2 * g1 * ((cos(w[k]) - (ahatx + (ahatx * bhatx * muhatx))) /
                                                (1 + (((ahatx + (ahatx * bhatx * muhatx))^2) - 2 * ((ahatx + (ahatx * bhatx * muhatx)) * cos(w[k]))))))
    
    if (f_X_wk <= 0) {
      return(Inf) # Penalize non-positive values to avoid log of zero or negative values
    }
    
    log_likelihood <- log_likelihood + (log(f_X_wk) + (Iprio[k] / f_X_wk))
  }
  return(log_likelihood / N) # Return normalized log-likelihood
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
# Calculating mu and sigma^2
#===============================================================================
# mu calculation based on estimated mu_eps
mu <- ((mu_eps + 2) / (mu_eps * (mu_eps + 1)))

# sigma^2 calculation with mu_eps
sigma2 <- (2 * (mu_eps + 2) / (mu_eps * (mu_eps + 1)) + 
             2 * (mu_eps + 3) / (mu_eps^2 * (mu_eps + 1)) - 
             ((mu_eps + 2) / (mu_eps * (mu_eps + 1)))^2)

#===============================================================================
# Prediction of the last ten data points using the whole data
#===============================================================================
# Truncate the dataset to the first 341 observations
x <- x[1:341]
N <- length(x)

#===============================================================================
# Residuals Calculation
#===============================================================================
# Initialize residuals vector e
e <- numeric(N)
e[1] <- 1  # Set the initial residual

# Compute residuals e[i] iteratively for i = 2 to N
for (i in 2:N) {
  e[i] <- max(0, x[i] - (a * x[i - 1]) - (a * b * x[i - 1] * e[i - 1]))
}

#===============================================================================
# Prediction of Future Values (xhat for h = 1 to h = 10)
#===============================================================================
# Initialize xhat for storing future predictions
xhat <- numeric(N + 10)  # Allocate space for N + 10 elements
xhat[N] <- x[N]          # Set the last observed value as the initial xhat

#===============================================================================
# Prediction for h = 1 (next step forecast)
#===============================================================================
xhat[N + 1] <- (a + (a * b * e[N])) * xhat[N] + mu  # First prediction

#===============================================================================
# Prediction for h = 2 to h = 10 (multi-step forecast)
#===============================================================================
# Loop to calculate future values for h = 2 to h = 10
for (h in 2:10) {
  sum_term <- 0  # Initialize sum term for the inner summation
  
  # Inner loop to calculate the sum term
  for (i in 0:(h - 2)) {
    sum_term <- sum_term + (a + a * b * mu)^i
  }
  
  # Calculate the predicted value for xhat[N + h]
  xhat[N + h] <- (a + a * b * mu)^(h - 1) * xhat[N + 1] + (a * b * sigma2 + mu) * sum_term
}

#===============================================================================
# Display the Predicted Values (xhat for the next 10 periods)
#===============================================================================
Xhat <- xhat[(N + 1):(N + 10)]  # Extract the last 10 predicted values

#===============================================================================
# End of Script
#===============================================================================
