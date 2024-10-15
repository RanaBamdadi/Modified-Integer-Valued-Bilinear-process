#===============================================================================
##################Residual analysis and model adequacy
#===============================================================================
################## PL-MINBL(1,0,1,1) model 
################## Poisson_lindly Modified Integer-Valued Bilinear
#===============================================================================
# Clear all variables in the environment
rm(list = ls(all=TRUE))

# Load necessary libraries
library(TSA)
library(stats)
library(rootSolve)
library(readxl)

#===============================================================================
# Read in Data
#===============================================================================
my_data <- read_excel("/Users/Documents/Data/RCI_offencebymonth.xlsm")

# Extract the data for time series analysis (specific rows and columns)
x <- as.numeric(my_data[3201, 4:354])  

# Additional vectors for periodogram and other calculations
Iprio <- c()
w <- c()

# Length of the data
N <- n <- length(x)

#===============================================================================
# Autocovariance Calculation
#===============================================================================
g0 <- acf(x, 50, type = "covariance", plot=FALSE)$acf[1]
g1 <- acf(x, 50, type = "covariance", plot=FALSE)$acf[2]
g2 <- acf(x, 50, type = "covariance", plot=FALSE)$acf[3]

#===============================================================================
# Nonlinear Function to Estimate Parameters 'a', 'b', and 'mu'
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
# Solve nonlinear system for parameter estimates using multiroot
ss <- multiroot(f = fn, start =c(0.284446,-0.159125,2.3623))

# Extract estimated parameters
aahatx=ss$root[1]
bbhatx=ss$root[2]
muhatx=ss$root[3]

#===============================================================================
# Periodogram Calculation
#===============================================================================
h <- (N-1)/2
Iprio <- periodogram(x, log='no', plot=FALSE)$spec[1:h]
w <- (2 * pi * (1:h)) / (N-1)

#===============================================================================
# Whittle Likelihood Function Definition
#===============================================================================
whittle_likelihood <- function(l) {
  ahatx <- l[1]
  bhatx <- l[2]
  muhatx <- (l[3]+2)/(l[3]*(l[3]+1))
  log_likelihood <- 0
  
  for (k in 1:h) {
    f_X_wk <- (1 / (2 * pi)) *(g0 + 2 * g1 * ((cos(w[k]) - (ahatx + (ahatx * bhatx * muhatx))) /
                                                (1 + (((ahatx + (ahatx * bhatx * muhatx))^2) - 2 * ((ahatx + (ahatx * bhatx * muhatx)) * cos(w[k]))))))
    
    if (f_X_wk <= 0) {
      return(Inf) # Penalize non-positive values to avoid log of zero or negative values
    }
    
    log_likelihood <- log_likelihood + (log(f_X_wk) + (Iprio[k] / f_X_wk))
  }
  return(log_likelihood / N)
}

#===============================================================================
# Parameter Estimation using Whittle Likelihood
#===============================================================================
initial_params <- c(aahatx,bbhatx,muhatx)

# Estimate parameters using nlm
est <- nlm(whittle_likelihood, initial_params)

# Extract estimated parameters
a <- est$estimate[1]
b <- est$estimate[2]
m <- est$estimate[3]

#===============================================================================
# Define parameters and initialize vectors for residual calculations
#===============================================================================

# Set the estimated parameters from the previous Whittle likelihood estimation
theta <- c(a, b, m)

# Define parameters and initialize vectors for residual calculations
n <- length(x)
E <- rep(NA, n)  # Expectations
V <- rep(NA, n)  # Variances
R <- c()         # Residuals
e <- c()         # Errors
lan <- c()       # Latent variables
p <- theta[3] / (1 + theta[3])  # Probability based on theta
u <- runif(n)    # Generate uniform random numbers

#===============================================================================
# Error Simulation based on Poisson-Lindley distribution
#===============================================================================

for (i in 1:n) {
  if (u[i] < p) {
    lan[i] <- rexp(1, theta[3])  # Exponential distribution with rate theta[3]
    e[i] <- rpois(1, lan[i])     # Poisson distribution with lambda from Exponential
  } else {
    lan[i] <- rexp(1, theta[3]) + rexp(1, theta[3])  # Sum of two Exponentials
    e[i] <- rpois(1, lan[i])     # Poisson distribution with combined lambda
  }
}

#===============================================================================
# Expectation and Variance Calculation
#===============================================================================

E[1] <- 1
V[1] <- 1  # Initialize for t = 1

for (t in 2:n) {
  # Expectation calculation
  E[t] <- (theta[1] + theta[1] * theta[2] * e[t - 1]) * x[t - 1] + 
    ((theta[3] + 2) / (theta[3] * (theta[3] + 1)))
  
  # Variance calculation
  V[t] <- theta[1]^2 * abs(theta[2]) * (1 - abs(theta[2])) * x[t - 1]^2 * abs(e[t - 1]) + 
    abs(theta[1]) * (1 - abs(theta[1])) * abs(x[t - 1]) + 
    abs(theta[1]) * (1 - abs(theta[1])) * theta[2]^2 * e[t - 1]^2 * abs(x[t - 1]) + 
    abs(theta[1]) * (1 - abs(theta[1])) * abs(theta[2]) * (1 - abs(theta[2])) * abs(e[t - 1]) * abs(x[t - 1]) + 
    2 * abs(theta[1]) * (1 - abs(theta[1])) * theta[2] * e[t - 1] * abs(x[t - 1]) + 
    2 * ((theta[3] + 2) / (theta[3] * (theta[3] + 1))) + 
    2 * ((theta[3] + 3) / ((theta[3]^2) * (theta[3] + 1))) - 
    (((theta[3] + 2) / (theta[3] * (theta[3] + 1)))^2)
}

#===============================================================================
# Standardize residuals
#===============================================================================

for (i in 1:n) {
  R[i] <- (x[i] - E[i]) / sqrt(V[i])  # Standardize residuals
}

#===============================================================================
# Plot Cumulative Periodogram of residuals and ACF
#===============================================================================

par(mfrow = c(1, 2))  # Set up 1x2 plotting layout

# Plot Cumulative Periodogram of residuals
cpgram(R, main = "")

# Plot ACF (Autocorrelation Function) of residuals
acf(R, main = "")

#===============================================================================
# Perform Ljung-Box test on residuals
#===============================================================================
Box.test(R, type = "Ljung-Box")
