#===============================================================================
################## Forecasting study
# The code in this file reproduces the Bootsrap forecasting results 
#===============================================================================
############################################
################## PL-MINBL(1,0,1,1) model 
################## Poisson_lindly Modified Integer-Valued Bilinear
############################################
# Clear all variables in the environment
rm(list = ls(all=TRUE))

# Load necessary libraries
library(TSA)
library(stats)
library(rootSolve)
library(readxl)
library(EnvStats)

#===============================================================================
# Initialization of vectors
#===============================================================================
ax <- c()    # To store estimates of parameter 'a'
bx <- c()    # To store estimates of parameter 'b'
mux <- c()   # To store estimates of 'mu'
xb <- c()    # Vector to store bootstrapped x values
e <- c()     # Residual vector

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
fn <- function(a) {
  c(
    F1 = (a[1] * a[2] * (((2 * ((a[3] + 2) / (a[3] * (a[3] + 1)))) + 
                            (2 * ((a[3] + 3) / (((a[3])^2) * (a[3] + 1))))) - 
                           (((a[3] + 2) / (a[3] * (a[3] + 1)))^2))) + 
      ((a[3] + 2) / (a[3] * (a[3] + 1))) - 
      mean(x) * (1 - (a[1] + a[1] * a[2] * ((a[3] + 2) / (a[3] * (a[3] + 1))))),
    
    F2 = (a[1] + a[1] * a[2] * ((a[3] + 2) / (a[3] * (a[3] + 1)))) - g2 / g1,
    
    F3 = (a[1] * a[2] * ((6 * ((a[3] + 4) / (((a[3])^3) * (a[3] + 1)))) + 
                           (12 * ((a[3] + 3) / (((a[3])^2) * (a[3] + 1)))) + 
                           (5 * (a[3] + 2) / (a[3] * (a[3] + 1))))) + 
      (a[1] * a[2] * mean(x) * (((2 * (a[3] + 2)) / (a[3] * (a[3] + 1))) + 
                                  ((2 * (a[3] + 3) / ((a[3]^2) * (a[3] + 1)))) - 
                                  (((a[3] + 2) / (a[3] * (a[3] + 1)))^2))) + 
      (-3 * a[1] * a[2] * ((a[3] + 2) / (a[3] * (a[3] + 1))) * 
         (((2 * (a[3] + 2)) / (a[3] * (a[3] + 1))) + 
            ((2 * (a[3] + 3) / ((a[3]^2) * (a[3] + 1)))) - 
            (((a[3] + 2) / (a[3] * (a[3] + 1)))^2))) - 
      (a[1] * a[2] * (((a[3] + 2) / (a[3] * (a[3] + 1)))^3)) - 
      (g1 - ((a[1] + a[1] * a[2] * ((a[3] + 2) / (a[3] * (a[3] + 1)))) * g0))
  )
}

# Solve nonlinear system for parameter estimates using multiroot
ss <- multiroot(f = fn, start = c(0.284446, -0.159125, 2.3623))

# Extract estimated parameters
aahatx <- ss$root[1]
bbhatx <- ss$root[2]
muhatx <- ss$root[3]

#===============================================================================
# Periodogram Calculation
#===============================================================================
h <- (n - 1) / 2
for (j in 1:h) {
  w[j] <- (2 * pi * j) / (n - 1)
  Iprio[j] <- periodogram(x, log = 'no', plot = FALSE)$spec[j]
}

#===============================================================================
# Whittle Likelihood Function Definition
#===============================================================================
whittle_likelihood <- function(l) {
  ahatx <- l[1]
  bhatx <- l[2]
  muhatx <- (l[3] + 2) / (l[3] * (l[3] + 1))
  log_likelihood <- 0
  
  # Calculate the log-likelihood
  for (k in 1:h) {
    f_X_wk <- (1 / (2 * pi)) * (g0 + 2 * g1 * 
                                  ((cos(w[k]) - (ahatx + (ahatx * bhatx * muhatx))) / 
                                     (1 + (((ahatx + (ahatx * bhatx * muhatx))^2) - 
                                             2 * ((ahatx + (ahatx * bhatx * muhatx)) * cos(w[k]))))))
    
    if (f_X_wk <= 0) {
      return(Inf)  # Penalize invalid values
    }
    
    log_likelihood <- log_likelihood + (log(f_X_wk) + (Iprio[k] / f_X_wk))
  }
  return(log_likelihood / N)
}

#===============================================================================
# Parameter Estimation using Whittle Likelihood
#===============================================================================
initial_params <- c(aahatx, bbhatx, muhatx)
est <- nlm(whittle_likelihood, initial_params)

# Extract estimated parameters from nlm
a <- est$estimate[1]
b <- est$estimate[2]
m <- est$estimate[3]

#===============================================================================
# Residual Calculation
#===============================================================================
e[1] <- 1  # Initial residual value

# Loop to calculate residuals for the series
for (i in 2:n) {
  e[i] <- round(max(0, (x[i] - a * x[i-1] - a * b * x[i-1] * e[i-1])))
}

#===============================================================================
# Bootstrap Simulation
#===============================================================================
B <- 1000  # Number of bootstrap samples
s <- 0
N <- n <- length(x) - 10  # Adjust length to remove future data points

# Bootstrapping loop
while (s <= B) {
  xb[1] <- 1  # Initial bootstrap value
  
  # Generate bootstrap samples
  for (i in 2:n) {
    xb[i] <- sign(a) * sign(xb[i-1]) * rbinom(1, abs(xb[i-1]), abs(a)) + 
      (sign(a) * sign(xb[i-1]) * rbinom(1, abs(xb[i-1]), abs(a))) * 
      (sign(b) * sign(e[i-1]) * rbinom(1, abs(e[i-1]), abs(b))) + e[i]
  }
#===============================================================================
# Nonlinear Function to Estimate Parameters 'a', 'b', and 'mu' for bootstrap samples
#===============================================================================
  ftotal<-function(a)  
    c(
      F1=(a[1]*a[2]*(((2*((a[3] + 2)/(a[3]*(a[3] + 1)))) + (2*((a[3] + 3)/(((a[3])^2)*(a[3] + 1))))) - (((a[3] + 2)/(a[3]*(a[3] + 1)))^2)))+((a[3] + 2)/(a[3]*(a[3] + 1)))-mean (xb)*(1 - (a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))),
      F2=(a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))-g2/g1,
      F3=(a[1]*a[2]*((6*((a[3] +4)/(((a[3])^3)*(a[3]+ 1)))) + (12*((a[3] +3)/(((a[3])^2)*(a[3]+ 1)))) + (5*(a[3]+2)/(a[3]*(a[3]+ 1)))))+
        (a[1]*a[2]*mean (xb)*(((2*(a[3]+2))/(a[3]*(a[3]+1)))+((2*(a[3]+3)/((a[3]^2)*(a[3]+1))))-(((a[3]+2)/(a[3]*(a[3]+1)))^2)))+
        (-3*a[1]*a[2]*((a[3] +2)/(a[3]*(a[3] +1)))*(((2*(a[3]+2))/(a[3]*(a[3]+1)))+((2*(a[3]+3)/((a[3]^2)*(a[3]+1))))-(((a[3]+2)/(a[3]*(a[3]+1)))^2)))-
        (a[1]*a[2]*(((a[3] + 2)/(a[3]*(a[3] + 1)))^3)) - (g1 - ((a[1] +a[1]* a[2]*((a[3] + 2)/(a[3]*(a[3] + 1))))*g0)) 
    )
  # Solve nonlinear system for parameter estimates using multiroot
  ss <- multiroot(ftotal, c(0.19419,0.0108033,2.27064), useFortran = FALSE)
  
  # Extract bootstrap samples estimated parameters
  aaahatx=ss$root[1]
  bbbhatx=ss$root[2]
  mmuhatx=ss$root[3]
  
  #===============================================================================
  # Periodogram Calculation
  #===============================================================================
  h=(n-1)/2
  for(j in 1:h){
    w[j]<-(2*pi*j)/(n-1)
    Iprio[j]<-periodogram(x,log='no',plot=FALSE)$spec[j]
  }
  
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
  
  # Extract estimated parameters from optimization
  ax[s]<-(est$estimate[1])
  bx[s]<-(est$estimate[2])
  mux[s]<-(est$estimate[3])
  
  s <- s + 1  # Increment bootstrap counter
}

#===============================================================================
# Calculate Average Parameter Estimates from Bootstrap
#===============================================================================
ab <- mean(ax)
bb <- mean(bx)
mb <- mean(mux)

#===============================================================================
# K-Step Ahead Prediction
#===============================================================================
# Initialize variables for the next 10 bootstrap replications
xbone <- c(); xbtwo <- c(); xbthree <- c(); xbfore <- c()
xbfive <- c(); xbsix <- c(); xbsev <- c(); xbeig <- c()
xbnine <- c(); xbten <- c()

# Calculate xb for the next 10 bootstrap periods using the recursive relation
xbone <- sign(ab) * sign(xb[n]) * rbinom(1, abs(xb[n]), abs(ab)) + 
  (sign(ab) * sign(xb[n]) * rbinom(1, abs(xb[n]), abs(ab))) * 
  (sign(bb) * sign(e[n]) * rbinom(1, e[n], abs(bb))) + e[n + 1]

xbtwo <- sign(ab) * sign(xbone) * rbinom(1, abs(xbone), abs(ab)) + 
  (sign(ab) * sign(xbone) * rbinom(1, abs(xbone), abs(ab))) * 
  (sign(bb) * sign(e[n + 1]) * rbinom(1, e[n + 1], abs(bb))) + e[n + 2]

xbthree <- sign(ab) * sign(xbtwo) * rbinom(1, abs(xbtwo), abs(ab)) + 
  (sign(ab) * sign(xbtwo) * rbinom(1, abs(xbtwo), abs(ab))) * 
  (sign(bb) * sign(e[n + 2]) * rbinom(1, e[n + 2], abs(bb))) + e[n + 3]

xbfore <- sign(ab) * sign(xbthree) * rbinom(1, abs(xbthree), abs(ab)) + 
  (sign(ab) * sign(xbthree) * rbinom(1, abs(xbthree), abs(ab))) * 
  (sign(bb) * sign(e[n + 3]) * rbinom(1, e[n + 3], abs(bb))) + e[n + 4]

xbfive <- sign(ab) * sign(xbfore) * rbinom(1, abs(xbfore), abs(ab)) + 
  (sign(ab) * sign(xbfore) * rbinom(1, abs(xbfore), abs(ab))) * 
  (sign(bb) * sign(e[n + 4]) * rbinom(1, e[n + 4], abs(bb))) + e[n + 5]

xbsix <- sign(ab) * sign(xbfive) * rbinom(1, abs(xbfive), abs(ab)) + 
  (sign(ab) * sign(xbfive) * rbinom(1, abs(xbfive), abs(ab))) * 
  (sign(bb) * sign(e[n + 5]) * rbinom(1, e[n + 5], abs(bb))) + e[n + 6]

xbsev <- sign(ab) * sign(xbsix) * rbinom(1, abs(xbsix), abs(ab)) + 
  (sign(ab) * sign(xbsix) * rbinom(1, abs(xbsix), abs(ab))) * 
  (sign(bb) * sign(e[n + 6]) * rbinom(1, e[n + 6], abs(bb))) + e[n + 7]

xbeig <- sign(ab) * sign(xbsev) * rbinom(1, abs(xbsev), abs(ab)) + 
  (sign(ab) * sign(xbsev) * rbinom(1, abs(xbsev), abs(ab))) * 
  (sign(bb) * sign(e[n + 7]) * rbinom(1, e[n + 7], abs(bb))) + e[n + 8]

xbnine <- sign(ab) * sign(xbeig) * rbinom(1, abs(xbeig), abs(ab)) + 
  (sign(ab) * sign(xbeig) * rbinom(1, abs(xbeig), abs(ab))) * 
  (sign(bb) * sign(e[n + 8]) * rbinom(1, e[n + 8], abs(bb))) + e[n + 9]

xbten <- sign(ab) * sign(xbnine) * rbinom(1, abs(xbnine), abs(ab)) + 
  (sign(ab) * sign(xbnine) * rbinom(1, abs(xbnine), abs(ab))) * 
  (sign(bb) * sign(e[n + 9]) * rbinom(1, e[n + 9], abs(bb))) + e[n + 10]

                                                         