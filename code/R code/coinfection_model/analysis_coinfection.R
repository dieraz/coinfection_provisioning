# Load required libraries for MCMC diagnostics, plotting, and data manipulation
library('fitR')         # For model fitting and diagnostics
library('MASS')         # For statistical tools (e.g., multivariate analysis)
library(tidyr)          # For tidy data reshaping
library(ggplot2)        # For plotting
library('coda')         # For MCMC diagnostics (trace plots, convergence)
library(lattice)        # For multi-panel plotting (e.g., density/trace plots)
library(reshape2)       # For reshaping data for plotting
library(plyr)           # For manipulating and combining data frames
library(stringr)        # For string operations
library(gridExtra)      # For arranging multiple plots in a grid

# Load required data and model definitions
load("Data_coinfection.RData")           # Contains prevalence data 
source('coinf_det_I9.R')        # Load deterministic co-infection model (fitmodel1)
source('coinf_det_L9.R')        # Load deterministic co-infection model (fitmodel2)
source('functions.R')           # Load custom functions (likelihoods, plotting, etc.)

# Load first MCMC chain results
load("COINF.ch1.RData")
my_mcmc.TdC2$acceptanceRate     # Print MCMC acceptance rate

# Process MCMC trace: convert to coda object, burn-in, thinning, summarise
trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 10000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 10)
theta1 <- colMeans(my_mcmc.TdC2$trace[,])  # Not using burn-in here

# Repeat for second chain
load("COINF.ch2.RData")
my_mcmc.TdC2$acceptanceRate
trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 10000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 10)
theta2 <- colMeans(my_mcmc.TdC2$trace[,])

# Repeat for third chain
load("COINF.ch3.RData")
my_mcmc.TdC2$acceptanceRate
trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 10000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 10)
theta3 <- colMeans(my_mcmc.TdC2$trace[,])

# Repeat for fourth chain
load("COINF.ch4.RData")
my_mcmc.TdC2$acceptanceRate
trace4 <- mcmc(my_mcmc.TdC2$trace)
trace.burn4 <- burnAndThin(trace4, burn = 10000)
trace.burn.thin4 <- burnAndThin(trace.burn4, thin = 10)
theta4 <- colMeans(my_mcmc.TdC2$trace[,])

# Combine thinned traces from all four MCMC chains into a coda object
trace.info <- mcmc.list(trace.burn.thin1, trace.burn.thin2, trace.burn.thin3, trace.burn.thin4)

# Visualize posterior distributions and convergence
densityplot(trace.info)       # Plot posterior densities
gelman.diag(trace.info)       # Gelman-Rubin statistic to check convergence

# Combine all traces into a single data frame for summarization
trace.combined <- ldply(trace.info)

# Compute posterior mean estimates for parameters across all chains
theta.bar <- colMeans(trace.combined[Coinf_det_w$thetaNames])

# Set initial state for simulation (host and parasite compartments)
init.state <- c(S = 80, I = 1, E = 0, L = 0, Ps = 1, Pi = 0)

# Compute log-likelihood for model 1 (prevE_HaddonW, e.g., environmental egg data)
log.like.theta.bar1 <- dTrajObs_DE6(Coinf_det_w, theta.bar, init.state, data = prevE_HaddonW, log = TRUE)
D.theta.bar1 <- -2 * log.like.theta.bar1

# Compute log-likelihood for model 2 (Hpoly_HaddonW, e.g., prevalence of another parasite stage)
log.like.theta.bar2 <- dTrajObs_DE4(Coinf_det_w2, theta.bar, init.state, data = Hpoly_HaddonW, log = TRUE)
D.theta.bar2 <- -2 * log.like.theta.bar2

# Compute effective number of parameters (p.D) using variance of -2*log posterior
p.D <- var(-2 * trace.combined$logDensity / 2)

# Calculate Deviance Information Criterion (DIC) for each model and total
DIC1 <- D.theta.bar1 + 2 * p.D
DIC2 <- D.theta.bar2 + 2 * p.D
DIC <- -2 * (log.like.theta.bar1 + log.like.theta.bar2) + 2 * p.D  # Joint model DIC

# Plot posterior predictive simulations versus observed data
prevplot <- plotFit_DE(Coinf_det_w, theta.bar, init.state, data = prevE_HaddonW, n.replicates = 20)
eggplot <- plotFit_DE(Coinf_det_w2, theta.bar, init.state, data = Hpoly_HaddonW, n.replicates = 20)

# Load Rmisc for computing confidence intervals
library(Rmisc)

# Compute 95% credible intervals for parameters
CI(trace.combined$delta, ci = 0.95)     
CI(trace.combined$gamma0, ci = 0.95)
CI(trace.combined$beta, ci = 0.95)
CI(trace.combined$v, ci = 0.95)
CI(trace.combined$gamma, ci = 0.95)
CI(trace.combined$sigma, ci = 0.95)





