# -----------------------------------------------------------
# Simulation and Parameter Estimation By Yule-Walker 
# -----------------------------------------------------------

# Clear the environment and load necessary libraries
rm(list = ls(all = TRUE))      # Clear all existing variables
library(itsmr)                 # For time series modeling
library(tseries)               # Time series analysis and related functions
library(Metrics)               # For performance evaluation (bias, MSE)
library(MLmetrics)             # For additional metrics like MSE

# ---------------------------------------------------------------------------------
# Initialize storage vectors for parameter estimates
# ---------------------------------------------------------------------------------
alpha_estimates <- c()         # Store estimates of alpha
beta_estimates <- c()          # Store estimates of beta
mu_estimates <- c()            # Store estimates of mu
Ahat <- c()                    # Temporary storage for A estimates
Bhat <- c()                    # Temporary storage for B estimates

# Set up simulation loop parameters
num_simulations <- 100         # Number of simulations to run
h <- 1                         # Simulation counter

# ---------------------------------------------------------------------------------
# Simulation loop: Iterate until num_simulations is reached
# ---------------------------------------------------------------------------------
while (h <= num_simulations) {
  
  # ---------------------------------------------------------------------------------
  # Step 1: Set parameters and generate initial data
  # ---------------------------------------------------------------------------------
  n <- 1001                     # Number of observations
  alpha <- -0.2                 # True alpha value
  beta <- -0.3                  # True beta value
  mu <- 7                       # True mu value
  
  x <- rep(NA, n)               # Initialize process vector
  e <- rpois(n, mu)             # Poisson-distributed errors (noise)
  y <- rep(NA, n)               # Initialize y process
  
  # Set initial value of the process
  y[1] <- rpois(1, abs(alpha))  # First value from Poisson
  
  # ---------------------------------------------------------------------------------
  # Step 2: Simulate the process
  # ---------------------------------------------------------------------------------
  for (i in 2:n) {
    y[i] <- sign(alpha) * sign(y[i - 1]) * rbinom(1, abs(y[i - 1]), abs(alpha)) +
      (sign(alpha) * sign(y[i - 1]) * rbinom(1, abs(y[i - 1]), abs(alpha))) * 
      (sign(beta) * sign(e[i - 1]) * rbinom(1, abs(e[i - 1]), abs(beta))) +
      e[i]
  }
  
  # Trim first value for analysis
  x <- y[2:n]
  e <- e[2:n]
  
  # ---------------------------------------------------------------------------------
  # Step 3: Compute autocovariances
  # ---------------------------------------------------------------------------------
  g0 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[1]   # Autocovariance at lag 0
  g1 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[2]   # Autocovariance at lag 1
  g2 <- acf(x, 50, type = "covariance", plot = FALSE)$acf[3]   # Autocovariance at lag 2
  
  # ---------------------------------------------------------------------------------
  # Step 4: Estimate parameters using Yule-Walker equations
  # ---------------------------------------------------------------------------------
  Ahat[h] <- g2 / g1                                    # Estimate for A
  Bhat[h] <- (g1 - (Ahat[h] * g0)) / (mean(x) + 1)      # Estimate for B
  mu_estimates[h] <- mean(x) * (1 - Ahat[h]) - Bhat[h]  # Estimate for mu
  alpha_estimates[h] <- Ahat[h] - Bhat[h]               # Estimate for alpha
  beta_estimates[h] <- Bhat[h] / (alpha_estimates[h] * mu_estimates[h])  # Estimate for beta
  
  # ---------------------------------------------------------------------------------
  # Step 5: Check if estimated parameters satisfy conditions
  # ---------------------------------------------------------------------------------
  if (-alpha_estimates[h] > 0 && -alpha_estimates[h] < 1 &&
      -beta_estimates[h] > 0 && -beta_estimates[h] < 1 &&
      mu_estimates[h] < mu + 0.5 && mu_estimates[h] > mu - 0.5) {
    
    # Increment simulation counter if conditions are satisfied
    h <- h + 1
  } else {
    # If conditions are not satisfied, continue to the next iteration
    print("Estimated parameters do not meet the required conditions.")
  }
}

# ---------------------------------------------------------------------------------
# Step 6: Calculate Bias and MSE for the parameter estimates
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