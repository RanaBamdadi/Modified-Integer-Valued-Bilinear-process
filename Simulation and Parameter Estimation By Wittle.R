# -----------------------------------------------------------
# Simulation and Parameter Estimation By Wittle
# -----------------------------------------------------------

# Clear the environment and load necessary libraries
rm(list = ls(all = TRUE))   # Clean workspace
library(TSA)                # For time series analysis
library(stats)              # Base statistics functions
library(Metrics)            # For bias and MSE calculations
library(MLmetrics)          # Additional metrics

# ---------------------------------------------------------------------------------
# Initialize parameters and storage vectors
# ---------------------------------------------------------------------------------
alpha_estimates <- c()      # Vector to store estimates of alpha
beta_estimates <- c()       # Vector to store estimates of beta
mu_estimates <- c()         # Vector to store estimates of mu

num_simulations <- 100      # Number of simulations to run
t <- 1                      # Counter for simulations

# ---------------------------------------------------------------------------------
# Simulation loop: Repeat the process 'num_simulations' times
# ---------------------------------------------------------------------------------
while (t <= num_simulations) {
  
  # Set parameters for the process
  N <- 1001                # Number of observations
  alpha <- -0.2            # True value of alpha
  beta <- -0.3             # True value of beta
  mu <- 7                  # True value of mu
  x <- rep(NA, N)          # Initialize time series vector
  
  # Generate Poisson-distributed errors
  epsilon <- rpois(N, mu)  # Poisson noise with mean 'mu'
  
  # Set the initial value of the time series
  x[1] <- 1                # Initial condition
  
  # ---------------------------------------------------------------------------------
  # Simulate the process based on the specified model
  # ---------------------------------------------------------------------------------
  for (i in 2:N) {
    x[i] <- sign(alpha) * sign(x[i - 1]) * rbinom(1, abs(x[i - 1]), abs(alpha)) +
      (sign(alpha) * sign(x[i - 1]) * sign(beta) * sign(epsilon[i - 1]) *
         rbinom(1, abs(x[i - 1]), abs(alpha))) * 
      rbinom(1, epsilon[i - 1], abs(beta)) + epsilon[i]
  }
  
  # Remove the first value for analysis
  x <- x[2:N]
  
  # ---------------------------------------------------------------------------------
  # Compute autocovariances using the acf function
  # ---------------------------------------------------------------------------------
  g0 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[1]
  g1 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[2]
  g2 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[3]
  
  # ---------------------------------------------------------------------------------
  # Periodogram calculation for spectral analysis
  # ---------------------------------------------------------------------------------
  h <- (N - 1) / 2
  w <- c()                 # Vector to store frequencies
  I_prio <- c()            # Vector to store periodogram values
  
  for (j in 1:h) {
    w[j] <- (2 * pi * j) / (N - 1)
    I_prio[j] <- periodogram(x, log = 'no', plot = FALSE)$spec[j]
  }
  
  # ---------------------------------------------------------------------------------
  # Define the objective function for parameter estimation
  # ---------------------------------------------------------------------------------
  objective_function <- function(params) {
    alpha_hat <- params[1]
    beta_hat <- params[2]
    mu_hat <- params[3]
    log_likelihood <- 0
    
    for (k in 1:h) {
      spectrum_estimate <- g0 + 2 * g1 * ((cos(w[k]) - (alpha_hat + alpha_hat * beta_hat * mu_hat)) /
                                            (1 + ((alpha_hat + alpha_hat * beta_hat * mu_hat)^2) - 
                                               2 * (alpha_hat + alpha_hat * beta_hat * mu_hat) * cos(w[k])))
      
      log_likelihood <- log_likelihood + (1 / N) * (log(spectrum_estimate) + I_prio[k] / spectrum_estimate)
    }
    
    return(log_likelihood)
  }
  
  # ---------------------------------------------------------------------------------
  # Estimate parameters using non-linear minimization
  # ---------------------------------------------------------------------------------
  initial_params <- c(alpha, beta, mu)   # Initial guess for parameters
  estimation <- nlm(objective_function, initial_params)
  
  # ---------------------------------------------------------------------------------
  # Check the estimated parameters meet the predefined conditions
  # ---------------------------------------------------------------------------------
  if (-estimation$estimate[1] > 0 && -estimation$estimate[1] < 1 &&
      -estimation$estimate[2] > 0 && -estimation$estimate[2] < 1 &&
      estimation$estimate[3] > mu - 0.5 && estimation$estimate[3] < mu + 0.5) {
    
    alpha_estimates[t] <- estimation$estimate[1]
    beta_estimates[t] <- estimation$estimate[2]
    mu_estimates[t] <- estimation$estimate[3]
    
    t <- t + 1  # Move to the next simulation
  } else {
    print("Estimated parameters do not meet the conditions.")
  }
}

# ---------------------------------------------------------------------------------
# Calculate Bias and MSE for parameter estimates
# ---------------------------------------------------------------------------------
cat("Bias and MSE for Alpha:\n")
cat("Bias: ", bias(alpha, alpha_estimates), "\n")
cat("MSE: ", MSE(alpha_estimates, alpha), "\n")

cat("Bias and MSE for Beta:\n")
cat("Bias: ", bias(beta, beta_estimates), "\n")
cat("MSE: ", MSE(beta_estimates, beta), "\n")

cat("Bias and MSE for Mu:\n")
cat("Bias: ", bias(mu, mu_estimates), "\n")
cat("MSE: ", MSE(mu_estimates, mu), "\n")