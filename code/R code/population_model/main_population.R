# Load required libraries for data manipulation, statistical modeling, and visualization
library('fitR')      # Package for fitting models to data
library('MASS')      # Contains functions for statistical methods, such as generating multivariate distributions
library(tidyr)       # Data wrangling library for reshaping data
library(ggplot2)     # Library for creating plots and data visualizations
library('coda')      # Library for MCMC diagnostics
library(deSolve)      # Solving differential equations (used for dynamic models)
library(lattice)     # Visualization library for creating trellis graphs

# Load custom functions and deterministic population dynamics model from external scripts
source('functions.R')  # Source custom R functions
source('popdyn_det.R')    # Source the population dynamics model script

# Load population data from an RData file
load("Data_population.RData")  # Contains population data 

# Fix the population data and remove incomplete cases (missing values)
POP_Haddon <- fix_POP_Haddon
POP_Haddon <- POP_Haddon[complete.cases(POP_Haddon), ]

# Set initial state and parameter values for the population dynamics model
init.state <- c(N = 50)  # Initial population size
theta <- c(b1 = 0.5, K = 50, mu = 1/24, w = 0.15)  # Initial parameter values

# Use the loaded population data
data1 <- POP_Haddon

# Define a function to compute the log-posterior probability of the model
posTdC <- function(theta) {
  my_fitmodel1 <- Popdyn_det       # Use the deterministic population dynamics model
  my_init.state <- init.state      # Set the initial state
  # Calculate the log-posterior for the given parameters using the population dynamics model
  return(dLogPosPOPDYN(fitmodel1 = my_fitmodel1,
                       theta = theta,
                       init.state = my_init.state,
                       data1 = data1))
}

# Set initial parameter values and bounds for the MCMC sampling
init.theta <- theta
lower <- c(b1 = 0.5, K = 25, mu = 1/52, w = 0)  # Lower bounds for parameters
upper <- c(b1 = 20, K = 100, mu = 1/8, w = 1)   # Upper bounds for parameters

# MCMC settings for the initial run
n.iterations <- 1000            # Number of iterations for MCMC
adapt.size.start <- 100         # Start of adaptive scaling of proposal distribution
adapt.size.cooling <- 0.999     # Cooling rate for the adaptation
adapt.shape.start <- 200        # Start of adaptive shaping of the proposal distribution

# Perform MCMC sampling using Metropolis-Hastings algorithm
my_mcmc.TdC <- mcmcMh(target = posTdC,           # Target log-posterior function
                      initTheta = theta,        # Initial parameter values
                      # proposal.sd = proposal.sd, # Proposal standard deviation (not used here)
                      limits = list(lower = lower, upper = upper), # Parameter bounds
                      nIterations = n.iterations,                 # Number of iterations
                      adaptSizeStart = adapt.size.start,         # Adaptive scaling start
                      adaptSizeCooling = adapt.size.cooling,     # Cooling rate for adaptation
                      adaptShapeStart = adapt.shape.start)       # Adaptive shaping start

# Update parameter estimates using the mean of the MCMC trace (excluding burn-in period)
theta <- colMeans(my_mcmc.TdC$trace[100:1000, 1:4])  # Estimate parameters based on MCMC results
covmat <- my_mcmc.TdC$covmat.empirical               # Obtain empirical covariance matrix

# Update MCMC settings for a more extensive sampling run
n.iterations <- 100000         # Increase the number of iterations
adapt.size.start <- 200        # Start of adaptive scaling
adapt.size.cooling <- 0.999    # Cooling rate remains the same
adapt.shape.start <- 200       # Start of adaptive shaping remains the same

# Perform a second, longer MCMC sampling run using updated parameter estimates and covariance matrix
my_mcmc.TdC2 <- mcmcMH(target = posTdC,               # Target log-posterior function
                       initTheta = theta,            # Updated initial parameter values
                       covmat = covmat,               # Use empirical covariance matrix from previous run
                       limits = list(lower = lower, upper = upper), # Parameter bounds
                       n.iterations = n.iterations,                  # Number of iterations
                       adaptSizeStart = adapt.size.start,          # Adaptive scaling start
                       adaptSizeCooling = adapt.size.cooling,      # Cooling rate for adaptation
                       adaptShapeStart = adapt.shape.start)        # Adaptive shaping start

# Save the results of the second MCMC sampling run to an RData file
save(my_mcmc.TdC2, file = 'POPDYN.1007.ch1.RData')
