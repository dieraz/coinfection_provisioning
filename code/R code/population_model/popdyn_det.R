# Define a function to simulate population dynamics based on given parameters, initial state, and times
simulate <- function(theta, init.state, times) {
  # Inner function defining the system of differential equations (ODEs)
  ode <- function(t, x, params) {
    
    # Extract the current population size
    N <- x[1]
    
    # Use the parameters and time to calculate the birth rate b
    with(as.list(params), {
      # The birth rate b is calculated using a sinusoidal function with a period of 52 (weeks)
      # This represents seasonal variation in the birth rate
      b <- abs(b1 * sin(2 * pi * (t / 52 - w))) + b1 * sin(2 * pi * (t / 52 - w))
      
      # Differential equation for the change in population size (dN)
      dN <- b * N * (K - N) / K - mu * N # Logistic growth with seasonal birth rate and constant mortality
      
      # Return the rate of change
      dx <- c(dN)
      list(dx)
    })
  }
  
  # Simulate the ODE using the lsoda solver, which handles stiff and non-stiff systems
  # Convert the output to a data frame
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  
  # Return the simulated trajectory
  return(traj)
}

# Define a function to generate a point observation based on the model's prediction
rPointObs <- function(model.point, theta) {
  # Generate a random observation from a Poisson distribution with the mean equal to the predicted population size
  obs.point <- rpois(n = 1, lambda = model.point[["N"]])
  
  # Return the observed data as a named vector
  return(c(obs = obs.point))
}

# Define a function to calculate the likelihood of observed data given the model prediction
dPointObs <- function(data.point, model.point, theta, log = FALSE) {
  # Compute the probability density of the observed data point based on the Poisson distribution
  # The mean of the distribution is given by the predicted population size
  return(dpois(x = data.point[["obs"]], lambda = model.point[["N"]], log = log))
}

# Define a function to calculate the prior distribution for the model parameters
dprior <- function(theta, log = FALSE) {
  # Calculate the log of uniform prior probabilities for each parameter
  log.prior.b1 <- dunif(theta[["b1"]], min = 0.5, max = 20, log = TRUE)   # Prior for b1
  log.prior.w <- dunif(theta[["w"]], min = 0, max = 1, log = TRUE)        # Prior for w
  log.prior.K <- dunif(theta[["K"]], min = 25, max = 100, log = TRUE)     # Prior for K
  log.prior.mu <- dunif(theta[["mu"]], min = 1/52, max = 1/8, log = TRUE) # Prior for mu
  
  # Sum the log-prior values
  log.sum = log.prior.b1 + log.prior.w + log.prior.K + log.prior.mu
  
  # Return the total prior probability in log form or as an exponentiated value
  return(ifelse(log, log.sum, exp(log.sum)))
}

# Set up the model using the fitmodel function
name <- "Population dynamics model"                  # Name of the model
state.names <- c("N")                                # State variable names (population size)
theta.names <- c("b1", "w", "K", "mu")               # Parameter names

# Create the population dynamics model object using the fitmodel function
Popdyn_det <- fitmodel(name, state.names, theta.names, simulate, rPointObs, dprior, dPointObs)
