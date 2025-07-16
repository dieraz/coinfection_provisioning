# Load required libraries for model fitting, statistical functions, and visualization
library('fitR')       # Package for fitting models to epidemiological data using MCMC
library('MASS')       # Functions for statistical modeling and multivariate distributions
library(tidyr)        # Data wrangling: reshaping and cleaning data
library(ggplot2)      # Visualization with grammar of graphics
library('coda')       # Diagnostics for MCMC output
library(deSolve)      # Solving differential equations (used for dynamic models)
library(lattice)      # Trellis plotting system

# Load external R scripts with custom model functions
source('functions.R')        # Contains likelihood and utility functions
source('coinf_det_I9.R')     # First deterministic co-infection model
source('coinf_det_L9.R')     # Second deterministic co-infection model

# Load empirical data
load("Data_coinfection.RData")           # Contains prevalence data 

# Set initial conditions for the compartmental model (S, I, E, L, Ps, Pi)
init.state <- c(S=70, I=1, E=0, L=0, Ps=1, Pi=0)

# Set initial parameter values for the co-infection model
theta <- c(
  delta = 1e-6,     # Interaction term or progression rate
  gamma0 = 1e-1,    # Recovery or transition rate
  beta = 1e-5,      # Transmission rate
  v = 1e-1,         # Re-infection or relapse rate
  gamma = 1e-1,     # Secondary recovery or transition rate
  sigma = 1e-1      # Possibly mortality or detection rate
)

# Assign empirical datasets to variables used in the likelihood function
data1 <- prevE_HaddonW       # Prevalence data for one pathogen
data2 <- Hpoly_HaddonW       # Data for the other co-infecting agent

# Define the log-posterior function combining two models and two datasets
posTdC <- function(theta) {
  my_fitmodel1 <- Coinf_det_w   # Co-infection model variant 1
  my_fitmodel2 <- Coinf_det_w2  # Co-infection model variant 2
  my_init.state <- init.state   # Initial state for model integration
  # Return the joint log-posterior using custom likelihood function
  return(dLogPos5(
    fitmodel1 = my_fitmodel1,
    fitmodel2 = my_fitmodel2,
    theta = theta,
    init.state = my_init.state,
    data1 = data1,
    data2 = data2
  ))
}

# Define prior bounds for each parameter (uniform priors)
init.theta <- theta
lower <- c(delta=1e-7, gamma0=1e-4, beta=1e-6, v=1e-3, gamma=1e-4, sigma=1e-3)
upper <- c(delta=1e-4, gamma0=1, beta=1e-1, v=1, gamma=1, sigma=1)

# MCMC settings for the short initial exploratory run
n.iterations <- 1000
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200

# Run MCMC sampling using Metropolis-Hastings algorithm
my_mcmc.TdC <- mcmcMh(
  target = posTdC,
  initTheta = theta,
  # proposal.sd = proposal.sd,   # (Not used — using adaptive proposals)
  limits = list(lower = lower, upper = upper),
  nIterations = n.iterations,
  adaptSizeStart = adapt.size.start,
  adaptSizeCooling = adapt.size.cooling,
  adaptShapeStart = adapt.shape.start
)

# Update parameter estimates using mean of MCMC samples (excluding burn-in)
theta <- colMeans(my_mcmc.TdC$trace[100:1000, 1:6])
covmat <- my_mcmc.TdC$covmat.empirical   # Empirical covariance of parameters

# Update MCMC settings for the long final sampling run
n.iterations <- 100000
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200


# Run final MCMC sampling using updated parameter values and empirical covariance
my_mcmc.TdC2 <- mcmcMh(
  target = posTdC,
  initTheta = theta,
  covmat = covmat,
  limits = list(lower = lower, upper = upper),
  nIterations = n.iterations,
  adaptSizeStart = adapt.size.start,
  adaptSizeCooling = adapt.size.cooling,
  adaptShapeStart = adapt.shape.start
)

# Save both initial and final MCMC objects to file
save(my_mcmc.TdC, my_mcmc.TdC2, file = 'COINF.ch1.RData')
