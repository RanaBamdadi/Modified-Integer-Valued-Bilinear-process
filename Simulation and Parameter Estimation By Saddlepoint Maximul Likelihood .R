# -----------------------------------------------------------
# Simulation and Parameter Estimation By Saddlepoint Maximul Likelihood 
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
    e<-rep(NA,n); e[1]=epsilon[1] # Initialize the first error
    d_epsilon_a<-rep(NA,n); d_epsilon_a[1]=0
    d_epsilon_b<-rep(NA,n); d_epsilon_b[1]=0
    d_epsilon_mu<-rep(NA,n); d_epsilon_mu[1]=0
    d_epsilon_a_a<-rep(NA,n); d_epsilon_a_a[1]=0
    d_epsilon_b_b<-rep(NA,n); d_epsilon_b_b[1]=0
    d_epsilon_mu_mu<-rep(NA,n); d_epsilon_mu_mu[1]=0
    d_epsilon_a_mu<-rep(NA,n); d_epsilon_a_b[1]=0
    d_epsilon_a_b<-rep(NA,n); d_epsilon_a_mu[1]=0
    d_epsilon_b_mu<-rep(NA,n); d_epsilon_b_mu[1]=0
    
    d_mu_a<-rep(NA,n); d_mu_a[1]=0
    d_mu_b<-rep(NA,n); d_mu_b[1]=0
    d_mu_mu<-rep(NA,n); d_mu_mu[1]=0
    d_mu_a_a<-rep(NA,n); d_mu_a_a[1]=0
    d_mu_b_b<-rep(NA,n); d_mu_b_b[1]=0
    d_mu_mu_mu<-rep(NA,n); d_mu_mu_mu[1]=0
    d_mu_a_b<-rep(NA,n); d_mu_a_b[1]=0
    d_mu_a_mu<-rep(NA,n); d_mu_a_mu[1]=0
    d_mu_b_mu<-rep(NA,n); d_mu_b_mu[1]=0
    
    d_sigma_a<-rep(NA,n); d_sigma_a[1]=0
    d_sigma_b<-rep(NA,n); d_sigma_b[1]=0
    d_sigma_mu<-rep(NA,n); d_sigma_mu[1]=0
    d_sigma_a_a<-rep(NA,n); d_sigma_a_a[1]=0
    d_sigma_b_b<-rep(NA,n); d_sigma_b_b[1]=0
    d_sigma_mu_mu<-rep(NA,n); d_sigma_mu_mu[1]=0
    d_sigma_a_b<-rep(NA,n); d_sigma_a_b[1]=0
    d_sigma_a_mu<-rep(NA,n); d_sigma_a_mu[1]=0
    d_sigma_b_mu<-rep(NA,n); d_sigma_b_mu[1]=0
    
    mu_t<-rep(NA,n); mu_t[1]=0
    sigma2_t<-rep(NA,n); sigma2_t[1]=0
    
    
    # Iterate over data points to calculate residuals and derivatives
    for(i in 2:n)
    {
      e[i]=round(max(0,(x[i]- params[1]*x[i-1]-params[1]*params[2]*x[i-1]*e[i-1]))) 
      d_epsilon_a[i]=(-(x[i-1])-(params[2]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a[i-1])) 
      d_epsilon_b[i]=(-( params[1]*x[i-1]*e[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b[i-1])) 
      d_epsilon_mu[i]=(-params[1]*params[2]*x[i-1]*d_epsilon_mu[i-1]) 
      d_epsilon_a_a[i]=(-(2*params[2]*x[i-1]*d_epsilon_a[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a_a[i-1])) 
      d_epsilon_b_b[i]=(-(2*params[1]*x[i-1]*d_epsilon_b[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b_b[i-1])) 
      d_epsilon_mu_mu[i]=(-params[1]*params[2]*x[i-1]*d_epsilon_mu_mu[i-1])
      d_epsilon_a_b[i]=((-x[i-1]*e[i-1])-(params[2]*x[i-1]*d_epsilon_b[i-1])-(params[1]*x[i-1]*d_epsilon_a[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a_b[i-1])) 
      d_epsilon_a_mu[i]=((-params[2]*x[i-1]*d_epsilon_mu[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_a_mu[i-1]))
      d_epsilon_b_mu[i]=((-params[1]*x[i-1]*d_epsilon_mu[i-1])-(params[1]*params[2]*x[i-1]*d_epsilon_b_mu[i-1]))
    }
    
    for (t in 2:n) 
    {
      d_mu_a[t] <- x[t-1] + params[2] * x[t-1] * e[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_a[t-1] 
      d_mu_b[t] <- params[1] * x[t-1] * e[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_b[t-1] 
      d_mu_mu[t] <- 1 + params[1] * params[2] * x[t-1] * d_epsilon_mu[t-1]
      d_mu_a_a[t] <- 2 * params[2] * x[t-1] * d_epsilon_a[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_a_a[t-1] 
      d_mu_b_b[t] <- 2 * params[1] * x[t-1] * d_epsilon_b[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_b_b[t-1] 
      d_mu_mu_mu[t] <- params[1] * params[2] * x[t-1] * d_epsilon_mu_mu[t-1] 
      d_mu_a_b [t]<- x[t-1] * e[t-1] + params[2] * x[t-1] * d_epsilon_b[t-1] + params[1] * x[t-1] * d_epsilon_a[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_a_b[t-1] 
      d_mu_a_mu[t] <- params[2] * x[t-1] * d_epsilon_mu[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_a_mu[t-1]
      d_mu_b_mu[t] <- params[1] * x[t-1] * d_epsilon_mu[t-1] + params[1] * params[2] * x[t-1] * d_epsilon_b_mu[t-1] 
    }
    
    for (t in 2:n) 
    {
      d_sigma_a[t] <- 
        2 * params[1]* abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * abs(e[t-1]) +
        (params[1]^2) * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a[t-1] +
        (sign(params[1]) - 2 * params[1]) * abs(x[t-1]) +
        (sign(params[1]) - 2 * params[1]) * (params[2])^2 * (e[t-1])^2 * abs(x[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (2 * e[t-1]) * d_epsilon_a[t-1] +
        (sign(params[1]) - 2 * params[1]) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * abs(e[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a[t-1] +
        2 * (sign(params[1]) - 2 * params[1]) * params[2] * abs(x[t-1]) * e[t-1] +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_a[t-1]
      
      d_sigma_b[t] <- 
        params[1]^2 * (sign(params[2]) - 2 * params[2]) * x[t-1]^2 * abs(e[t-1]) +
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_b[t-1] +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * e[t-1]^2 * abs(x[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (2 * e[t-1]) * d_epsilon_b[t-1] +
        abs(params[1]) * (1 - abs(params[1])) * (sign(params[2]) - 2 * params[2]) * abs(e[t-1]) * abs(x[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_b[t-1] +
        2 * abs(params[1]) * (1 - abs(params[1])) * e[t-1] * abs(x[t-1]) +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_b[t-1]
      
      d_sigma_mu[t] <- 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_mu[t-1] +
        abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (2 * e[t-1]) * d_epsilon_mu[t-1] +
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_mu[t-1] +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_mu[t-1]+1
      
      d_sigma_a_a[t] <-
        2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * abs(e[t-1]) +
        4 * params[1] * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a[t-1] +
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a_a[t-1] -
        2 * abs(x[t-1]) -
        2 * params[2]^2 * e[t-1]^2 * abs(x[t-1]) + 
        4 * (sign(params[1]) - 2 * params[1]) * params[2]^2 * abs(x[t-1]) * (e[t-1]) * d_epsilon_a[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * ( d_epsilon_a[t-1]*d_epsilon_a[t-1] + e[t-1] * d_epsilon_a_a[t-1]) - 
        2 * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) *abs(e[t-1]) + 
        2*(sign(params[1]) - 2 * params[1]) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a_a[t-1] - 
        4 * params[2] * abs(x[t-1]) * e[t-1] + 
        4 * (sign(params[1]) - 2 * params[1]) * params[2] * abs(x[t-1]) * d_epsilon_a[t-1] +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_a_a[t-1]
      
      d_sigma_b_b[t] <- 
        - 2 * params[1]^2 * x[t-1]^2 * abs(e[t-1])  + 
        2*params[1]^2 * (sign(params[2]) - 2 * params[2]) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_b[t-1] + 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_b_b[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * e[t-1]^2 * abs(x[t-1]) + 
        8 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * e[t-1] * d_epsilon_b[t-1] + 
        2 *abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * ( d_epsilon_b[t-1]*d_epsilon_b[t-1] +  e[t-1] * d_epsilon_b_b[t-1]) - 
        2 * abs(params[1]) * (1 - abs(params[1])) * abs(e[t-1]) * abs(x[t-1]) + 
        2*abs(params[1]) * (1 - abs(params[1])) * (sign(params[2]) - 2 * params[2]) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_b[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2]))*abs(x[t-1]) * sign(e[t-1]) * d_epsilon_b_b[t-1] + 
        4 * abs(params[1]) * (1 - abs(params[1])) * abs(x[t-1]) * d_epsilon_b[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_b_b[t-1]
      
      
      d_sigma_mu_mu[t] <- 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_mu_mu[t-1] + 
        2*abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (d_epsilon_mu[t-1]*d_epsilon_mu[t-1] + e[t-1] * d_epsilon_mu_mu[t-1]) + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_mu_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_mu_mu[t-1]
      
      d_sigma_a_b[t] <- 
        2 * params[1] * (sign(params[2]) - 2 * params[2]) * x[t-1]^2 * abs(e[t-1]) + 
        2 * params[1] * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_b[t-1] + 
        params[1]^2 * (sign(params[2]) - 2 * params[2]) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a[t-1] + 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a_b[t-1] + 
        2 * (sign(params[1]) - 2 * params[1]) * params[2] * e[t-1]^2 * abs(x[t-1]) + 
        2 * (sign(params[1]) - 2 * params[1]) * params[2]^2 * e[t-1] * abs(x[t-1]) * d_epsilon_b[t-1] + 
        4 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * e[t-1] * d_epsilon_a[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * ( d_epsilon_a[t-1]*d_epsilon_b[t-1] + e[t-1] * d_epsilon_a_b[t-1]) + 
        (sign(params[1]) - 2 * params[1]) * (sign(params[2]) - 2 * params[2]) * abs(x[t-1]) * abs(e[t-1]) + 
        (sign(params[1]) - 2 * params[1]) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_b[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * (sign(params[2]) - 2 * params[2]) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a_b[t-1] + 
        2 * (sign(params[1]) - 2 * params[1]) * abs(x[t-1]) * e[t-1] + 
        2 * (sign(params[1]) - 2 * params[1]) * params[2] * abs(x[t-1]) * d_epsilon_b[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * abs(x[t-1]) * d_epsilon_a[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_a_b [t-1]
      
      d_sigma_a_mu[t] <- 
        2 * params[1] * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_mu[t-1] + 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_a_mu[t-1] + 
        2 *(sign(params[1]) - 2 * params[1]) * params[2]^2 * abs(x[t-1]) *  e[t-1] * d_epsilon_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (d_epsilon_mu[t-1]*d_epsilon_a[t-1] + e[t-1] * d_epsilon_a_mu[t-1]) + 
        (sign(params[1]) - 2 * params[1]) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_mu[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_a_mu[t-1] + 
        2 * (sign(params[1]) - 2 * params[1]) * params[2] * abs(x[t-1]) * d_epsilon_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_a_mu[t-1]  
      
      d_sigma_b_mu[t] <- 
        params[1]^2 * (sign(params[2]) - 2 * params[2]) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_mu[t-1] + 
        params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * sign(e[t-1]) * d_epsilon_b_mu[t-1] + 
        4 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1])  * e[t-1]* d_epsilon_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * abs(x[t-1]) * (d_epsilon_mu[t-1]*d_epsilon_b[t-1] + e[t-1] * d_epsilon_b_mu[t-1]) + 
        abs(params[1]) * (1 - abs(params[1])) * (sign(params[2]) - 2 * params[2]) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_mu[t-1] + 
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(x[t-1]) * sign(e[t-1]) * d_epsilon_b_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * abs(x[t-1]) * d_epsilon_mu[t-1] + 
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * abs(x[t-1]) * d_epsilon_b_mu[t-1]
    }
    
    for (t in 2:n) {
      mu_t[t] <- (params[1] + params[1] * params[2] * e[t-1]) * x[t-1] + params[3]
      sigma2_t[t] <- params[1]^2 * abs(params[2]) * (1 - abs(params[2])) * x[t-1]^2 * abs(e[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * abs(x[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * params[2]^2 * e[t-1]^2 * abs(x[t-1]) +
        abs(params[1]) * (1 - abs(params[1])) * abs(params[2]) * (1 - abs(params[2])) * abs(e[t-1]) * abs(x[t-1]) +
        2 * abs(params[1]) * (1 - abs(params[1])) * params[2] * e[t-1] * abs(x[t-1]) +
        params[3]
    }
    
    i=2:n
    L_a <- sum(((x[i] - mu_t[i]) / (sigma2_t[i])) * d_mu_a[i]) + sum((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) * d_sigma_a[i])
    L_b <- sum(((x[i] - mu_t[i]) / (sigma2_t[i])) * d_mu_b[i] )+ sum( (((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) * d_sigma_b[i])
    L_mu <- sum(((x[i] - mu_t[i]) / (sigma2_t[i])) * d_mu_mu[i]) + sum((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) * d_sigma_mu[i]) 
    L_a_a <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_a_a[i]) - 2 * sum(((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_a[i] * d_sigma_a[i] )+ sum( ((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_a_a[i])
    L_b_b <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_b_b[i]) - 2 *sum( ((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_b[i] * d_sigma_b[i]) + sum(((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_b_b[i])
    L_mu_mu <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_mu_mu[i] )- 2 * sum(((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_mu[i] * d_sigma_mu[i] )+ sum(((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_mu_mu[i])
    L_a_b <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_a_b[i] )- 2 * sum(((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_a[i] * d_sigma_b[i] )+ sum(((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_a_b[i])
    L_a_mu <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_a_mu[i]) - 2 *sum (((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_a[i] * d_sigma_mu[i] )+ sum(((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_a_mu[i])
    L_b_mu <- sum((((x[i] - mu_t[i])-1 )/ (sigma2_t[i])) * d_mu_b_mu[i]) - 2 * sum(((x[i] - mu_t[i]) / (sigma2_t[i]^2)) * d_mu_b[i] * d_sigma_mu[i]) + sum(((((x[i] - mu_t[i])^2 - sigma2_t[i]) / (2 * sigma2_t[i]^2)) - (((x[i] - mu_t[i])^2) / (sigma2_t[i]^3)) + (1 / (2 * sigma2_t[i]^2))) * d_sigma_b_mu[i])
    
    
    G=c(L_a,L_b,L_mu)
    H=matrix(c(L_a_a,L_a_b,L_a_mu,L_a_b,L_b_b,L_b_mu,L_a_mu,L_b_mu,L_mu_mu),3,3)
    
    if (any(!is.finite(H))) 
    {
      print("Error: H matrix contains infinite values.")
      break
    }
    else{
      params_new <- c(params + (ginv(-H) %*% G))
      if (abs(params_new[1] - params[1]) < 10^(-3) && abs(params_new[2] - params[2]) < 10^(-3) && abs(params_new[3] - params[3]) < 10^(-3)) 
      {
        if ( -params[1] > 0 && -params[1] < 1 && -params[2] > 0  && -params[2] < 1 && params[3] > mu-0.5 && params[3] < mu+0.5 && -params_new[1] > 0  && -params_new[1] < 1 && -params_new[2] > 0  && -params_new[2] < 1  && params_new[3] > mu-0.5 && params_new[3] < mu+0.5) 
        {
          ahat[num_iterations] <- params_new[1]
          bhat[num_iterations] <- params_new[2]
          muhat[num_iterations] <- params_new[3]
          num_iterations<-num_iterations+1
          break
        }
        else{
          print( "Conditions related to params_new values not satisfied.")}
        break 
      } 
      else
      {
        params<-params_new
        j<-j+1
      }
    }
  }
  
}

# Calculate performance metrics
cat("Bias and MSE for parameter 'a': ", bias(a, ahat), MSE(ahat, a), "\n")
cat("Bias and MSE for parameter 'b': ", bias(b, bhat), MSE(bhat, b), "\n")
cat("Bias and MSE for parameter 'mu': ", bias(m, muhat), MSE(muhat, m), "\n")


