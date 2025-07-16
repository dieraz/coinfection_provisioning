# Load required libraries for data manipulation, plotting, and statistical analysis
library('fitR')      # For fitting models
library('MASS')      # For statistical functions
library(tidyr)       # For data tidying
library(ggplot2)     # For data visualization
library('coda')      # For Markov Chain Monte Carlo (MCMC) diagnostics
library(lattice)     # For plotting MCMC trace plots
library(reshape2)    # For data reshaping
library(plyr)        # For data manipulation functions
library(stringr)     # For string manipulation
library(gridExtra)   # For arranging multiple grid-based plots

# Load required data and model definitions
source('popdyn_det.R')             # Load population dynamics model
load("Data_population.RData")  # Contains population data 
load("POPDYN.1007.ch1.RData")      # Load MCMC results for the first chain

# Check the acceptance rate for the first MCMC chain
my_mcmc.TdC2$acceptance.rate

# Convert the trace of the MCMC results to an MCMC object for diagnostics
trace1 <- mcmc(my_mcmc.TdC2$trace)

# Plot the trace of the MCMC samples to assess convergence
xyplot(trace1)

# Perform burn-in by discarding the first 5000 samples to remove initial bias
trace.burn1 <- burnAndThin(trace1, burn = 5000)

# Thin the samples to reduce autocorrelation by keeping every 100th sample
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 100)

# Calculate the mean parameter estimates after burn-in and thinning
theta1 <- colMeans(trace.burn.thin1)

# Repeat the above steps for the second MCMC chain
load("POPDYN.1007.ch2.RData")
my_mcmc.TdC2$acceptance.rate
trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 100)
theta2 <- colMeans(trace.burn.thin2)

# Repeat the above steps for the third MCMC chain
load("POPDYN.1007.ch3.RData")
my_mcmc.TdC2$acceptance.rate
trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 100)
theta3 <- colMeans(trace.burn.thin3)

# Repeat the above steps for the fourth MCMC chain
load("POPDYN.1007.ch4.RData")
my_mcmc.TdC2$acceptance.rate
trace4 <- mcmc(my_mcmc.TdC2$trace)
trace.burn4 <- burnAndThin(trace4, burn = 5000)
trace.burn.thin4 <- burnAndThin(trace.burn4, thin = 100)
theta4 <- colMeans(trace.burn.thin4)

# Combine the thinned traces from all four MCMC chains for diagnostics
trace.info <- mcmc.list(trace.burn.thin1, trace.burn.thin2, trace.burn.thin3, trace.burn.thin4)

# Perform Gelman-Rubin diagnostic to assess convergence across the chains
gelman.diag(trace.info)

# Combine the traces into a single data frame
trace.combined <- ldply(trace.info)

# Calculate the mean parameter estimates across all combined chains
theta.bar <- colMeans(trace.combined[Popdyn_det$theta.names])

# Load custom functions for data analysis
source('my_functions.R')

# Preprocess the population data by removing rows with missing values
POP_Haddon <- fix_POP_Haddon
POP_Haddon <- POP_Haddon[complete.cases(POP_Haddon), ]

# Set the initial state for the population dynamics model
init.state <- c(N = 50)

# Calculate the log-likelihood of the data given the parameter estimates
log.like.theta.bar <- dTrajObs_POPDYN(Popdyn_det, theta.bar, init.state, data = POP_Haddon, log = TRUE)

# Calculate the deviance of the model (D_theta_bar)
D.theta.bar <- -2 * log.like.theta.bar

# Calculate the effective number of parameters (p.D) based on the variance of the log density
p.D <- var(-2 * trace.combined$log.density) / 2

# Calculate the Deviance Information Criterion (DIC) for model comparison
DIC <- D.theta.bar + 2 * p.D

# Plot the model fit against the data using the estimated parameters
plotFit_DE(Popdyn_det, theta.bar, init.state, data = POP_Haddon, n.replicates = 20)

# Load the Rmisc library for calculating confidence intervals
library(Rmisc)

# Calculate 95% confidence intervals for the parameter estimates from the combined MCMC samples
CI(trace.combined$b1, ci = 0.95)
CI(trace.combined$w, ci = 0.95)
CI(trace.combined$K, ci = 0.95)
CI(trace.combined$mu, ci = 0.95)
