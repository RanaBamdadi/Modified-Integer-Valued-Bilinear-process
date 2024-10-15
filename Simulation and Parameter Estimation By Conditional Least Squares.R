# -----------------------------------------------------------
# Simulation and Parameter Estimation By Conditional Least Squares
# -----------------------------------------------------------

# Clear environment and load necessary libraries
rm(list = ls(all=TRUE))  # Remove all objects from the environment
library(MASS)            # Matrix algebra
library(Metrics)         # Performance metrics
library(MLmetrics)       # Additional performance metrics

# Initialize vectors to store parameter estimates for 'a', 'b', and 'mu'
ahat <- c()  # Estimates for parameter 'a'
bhat <- c()  # Estimates for parameter 'b'
muhat <- c()  # Estimates for parameter 'mu'

# Set parameters for the simulation
a <- -0.2  # True value for parameter 'a'
b <- -0.3  # True value for parameter 'b'
m <- 7     # True value for mean 'mu'
N <- 1001  # Number of data points
num_iterations <- 1  # Counter for number of iterations

# Main simulation loop (repeat 100 iterations)
while(num_iterations <= 100) {
  
  # Generate Poisson-distributed errors and initialize process
  epsilon <- rpois(N, m)   # Error term following Poisson distribution
  h <- rep(NA, N)          # To store generated process
  h[1] <- rpois(1, abs(a)) # Initialize first value of the process
  
  # Simulate the process
  for (i in 2:N) {
    h[i] <- sign(a) * sign(h[i-1]) * rbinom(1, abs(h[i-1]), abs(a)) +
      (sign(a) * sign(h[i-1]) * rbinom(1, abs(h[i-1]), abs(a))) *
      (sign(b) * sign(epsilon[i-1]) * rbinom(1, abs(epsilon[i-1]), abs(b))) +
      epsilon[i]
  }
  
  # Remove the first value for further analysis
  x <- h[2:N]
  n <- length(x)
  
  # Initialize parameters for optimization
  params <- c(a, b, m)
  j <- 1  # Optimization loop counter
  
  # Optimization loop (max 20 iterations)
  while(j <= 20) {
    # Initialize residuals and derivatives
    e <- rep(NA, n); e[1] <- epsilon[1]
    d_epsilon_a <- rep(NA, n); d_epsilon_a[1] = 0
    d_epsilon_b <- rep(NA, n); d_epsilon_b[1] = 0
    d_epsilon_mu <- rep(NA, n); d_epsilon_mu[1] = 0
    d_epsilon_a_a <- rep(NA, n); d_epsilon_a_a[1] = 0
    d_epsilon_b_b <- rep(NA, n); d_epsilon_b_b[1] = 0
    d_epsilon_mu_mu <- rep(NA, n); d_epsilon_mu_mu[1] = 0
    d_epsilon_a_b <- rep(NA, n); d_epsilon_a_b[1] = 0
    d_epsilon_a_mu <- rep(NA, n); d_epsilon_a_mu[1] = 0
    d_epsilon_b_mu <- rep(NA, n); d_epsilon_b_mu[1] = 0
    
    # Compute residuals and derivatives
    for (i in 2:n) {
      e[i] <- round(max(0, (x[i] - params[1] * x[i-1] - params[1] * params[2] * x[i-1] * e[i-1])))
      d_epsilon_a[i] <- -(x[i-1] + params[2] * x[i-1] * e[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_a[i-1])
      d_epsilon_b[i] <- -(params[1] * x[i-1] * e[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_b[i-1])
      d_epsilon_mu[i] <- -params[1] * params[2] * x[i-1] * d_epsilon_mu[i-1]
      d_epsilon_a_a[i] <- -(2 * params[2] * x[i-1] * d_epsilon_a[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_a_a[i-1])
      d_epsilon_b_b[i] <- -(2 * params[1] * x[i-1] * d_epsilon_b[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_b_b[i-1])
      d_epsilon_mu_mu[i] <- -params[1] * params[2] * x[i-1] * d_epsilon_mu_mu[i-1]
      d_epsilon_a_b[i] <- -(x[i-1] * e[i-1] + params[2] * x[i-1] * d_epsilon_b[i-1] + params[1] * x[i-1] * d_epsilon_a[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_a_b[i-1])
      d_epsilon_a_mu[i] <- -(params[2] * x[i-1] * d_epsilon_mu[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_a_mu[i-1])
      d_epsilon_b_mu[i] <- -(params[1] * x[i-1] * d_epsilon_mu[i-1] + params[1] * params[2] * x[i-1] * d_epsilon_b_mu[i-1])
    }
    
    # Compute gradients and Hessians
    dqa=2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-x[i-1]-(params[2]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a[i-1])))
    dqb=2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-( params[1]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b[i-1])))
    dqm=2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-1-params[1]*params[2]*x[i-1]*d_epsilon_mu[i-1]))
    d2qa2=2*sum((-x[i-1]-(params[2]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a[i-1]))^2)+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-2*params[2]*x[i-1]*d_epsilon_a[i-1]-params[1]*params[2]*x[i-1]*d_epsilon_a_a[i-1]))
    d2qb2=2*sum((-( params[1]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b[i-1]))^2)+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-2*params[1]*x[i-1]*d_epsilon_b[i-1]-params[1]*params[2]*x[i-1]*d_epsilon_b_b[i-1]))
    d2qm2=2*sum((-1-params[1]*params[2]*x[i-1]*d_epsilon_mu[i-1])^2)+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-params[1]*params[2]*x[i-1]*d_epsilon_mu_mu[i-1]))
    d2qab=2*sum((-x[i-1]-(params[2]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a[i-1]))*(-( params[1]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b[i-1])))+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-(x[i-1]*e[i-1])-(params[2]*x[i-1]*d_epsilon_b[i-1])-(params[1]*x[i-1]*d_epsilon_a[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a_b[i-1])))
    d2qam=2*sum((-x[i-1]-(params[2]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a[i-1]))*(-1-params[1]*params[2]*x[i-1]*d_epsilon_mu[i-1]))+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-params[2]*x[i-1]*d_epsilon_mu[i-1]-params[1]*params[2]*x[i-1]*d_epsilon_a_mu[i-1]))
    d2qbm=2*sum((-( params[1]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b[i-1]))*(-1-params[1]*params[2]*x[i-1]*d_epsilon_mu[i-1]))+
      2*sum((x[i]-params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]-(params[3]))*(-params[1]*x[i-1]*d_epsilon_mu[i-1]-params[1]*params[2]*x[i-1]*d_epsilon_b_mu[i-1]))
    
    G=c(dqa,dqb,dqm)
    H=matrix(c(d2qa2,d2qab,d2qam,d2qab,d2qb2,d2qbm,d2qam,d2qbm,d2qm2),3,3)
    if (any(!is.finite(H))) 
    {
      print("Error: H matrix contains infinite or missing values.")
      break
    }
    else {
      
      params_new=c(params+((ginv(-H))%*%G))
      
      if(abs(params_new[1]-(params[1]))<10^(-3)&& abs(params_new[2]-(params[2]))<10^(-3)&& abs(params_new[3]-(params[3]))<10^(-3)) 
      {
        if ( -params[1] > 0 && -params[1] < 1 && -params[2] > 0  && -params[2] < 1 && params[3] > m-0.5 && params[3] < m+0.5 && -params_new[1] > 0  && -params_new[1] < 1 && -params_new[2] > 0  && -params_new[2] < 1  && params_new[3] > m-0.5 && params_new[3] < m+0.5) 
        {
          
          ahat[num_iterations]<-params_new[1]
          bhat[num_iterations]<-params_new[2]
          muhat[num_iterations]<-params_new[3]
          num_iterations<-num_iterations+1
          break
        }
        else{
          print( "Conditions related to params_new values not satisfied.")
        }
        break 
      } else{
        params<-params_new
        j<-j+1}
    }
  }
}
  

# Calculate performance metrics
cat("Bias and MSE for parameter 'a': ", bias(a, ahat), MSE(ahat, a), "\n")
cat("Bias and MSE for parameter 'b': ", bias(b, bhat), MSE(bhat, b), "\n")
cat("Bias and MSE for parameter 'mu': ", bias(m, muhat), MSE(muhat, m), "\n")
