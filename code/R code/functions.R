plotFit_DE <- function (fitmodel, theta, init.state, data, n.replicates = 1, 
          summary = TRUE, alpha = min(1, 10/n.replicates), all.vars = FALSE, 
          non.extinct = NULL, observation = TRUE, plot = TRUE){
  times <- rep(0:tail(data$time,n=1))
  if (n.replicates > 1) {
    cat("Simulate ", n.replicates, " replicate(s)\n")
  }
  traj <- simulateModelReplicates(fitmodel = fitmodel, theta = theta, 
                                  initState = init.state, times = times, n = n.replicates, 
                                  observation = observation)
  if (all.vars) {
    state.names <- NULL
  }
  else {
    state.names <- grep("obs", names(traj), value = TRUE)
  }
  p <- plotTraj(traj = traj, stateNames = state.names, data = data, 
                summary = summary, alpha = alpha, nonExtinct = non.extinct, 
                plot = FALSE)
  if (plot) {
    print(p)
  }
  else {
    return(list(traj = traj, plot = p))
  }
}

computeDistanceABC_DE <- function (sum.stats, distanceABC, fitmodel, theta, init.state,data){
  times <- rep(0:tail(data$time,n=1))
  
  model.obs <- rTrajObs(fitmodel = fitmodel, theta = theta, 
                        init.state = init.state, times = times)
  dist.ABC <- distanceABC(sum.stats = sum.stats, data.obs = data, 
                          model.obs = model.obs)
  return(dist.ABC)
}

dLogPos <- function(fitmodel, theta, init.state, data) {
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_DE(fitmodel, theta, init.state, data, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dTrajObs_DE <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["Eins"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}

dLogPos2 <- function(fitmodel, theta, init.state, data) {
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_DE2(fitmodel, theta, init.state, data, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dTrajObs_DE2 <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["prev"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}


dLogPos3 <- function(fitmodel, theta, init.state, data) {
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_DE3(fitmodel, theta, init.state, data, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dTrajObs_DE3 <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["Num"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}

dLogPos4 <- function(fitmodel, theta, init.state, data) {
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_DE4(fitmodel, theta, init.state, data, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dTrajObs_DE4 <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["Ldet"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}


dLogPos5 <- function(fitmodel1,fitmodel2, theta, init.state, data1, data2) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_DE6(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_DE4(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
  
}

dLogPos6 <- function(fitmodel, theta, init.state, data) {
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_DE6(fitmodel, theta, init.state, data, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dTrajObs_DE6 <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["I"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}



dLogPos.mult <- function(theta) {
  
  return(dLogPos(fitmodel = fitmodel,
                 theta = theta,
                 init.state = init.state,
                 data = data))
  
}

dLogPos5 <- function(fitmodel1,fitmodel2, theta, init.state, data1, data2) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_DE6(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_DE4(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
  
}

dTrajObs_DE_hf <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["prev"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}

dTrajObs_DE_hm <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["prev"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}



dTrajObs_DE_hsa <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["prev"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}

dLogPosherpes <- function(fitmodel1,fitmodel2, theta, init.state, data1, data2) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_DE_hf(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_DE_hm(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
  
}

dLogPosherpes2 <- function(fitmodel1,fitmodel2, fitmodel3, theta, init.state, data1, data2, data3) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_DE_hf(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_DE_hm(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior2 <- log.posterior1 + log.likelihood2
  
  log.likelihood3 <- dTrajObs_DE_hsa(fitmodel3, theta, init.state, data3, log = TRUE)
  log.posterior <- log.posterior2 + log.likelihood3
  
  return(log.posterior)
  
}

particleFilter_herpes <- function(fitmodel1,fitmodel2, fitmodel3, theta, init.state, data1, data2, data3, n.particles) {
  
  ## Initialisation of the algorithm
  
  # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
  margLogLike1 <- 0
  
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  
  # Initialise time variable
  current.time <- 0
  
  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data1))){
    
    # Extract next data1 point (must be a vector)
    data1.point <- unlist(data1[i, ])
    next.time <- data1.point["time"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    
    ## Loop over particles: propagate and weight
    for(p in 1:n.particles){
      
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel1$simulate`
      traj <- fitmodel1$simulate(theta = theta,
                                 init.state = current.state.particle,
                                 times = c(current.time,next.time))
      
      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector
      model.point <- unlist(traj[2,fitmodel1$state.names])
      
      # Weight the particle with the likelihood of the observed 
      # data1 point using the function `fitmodel1$dPointObs`
      weight.particles[p] <-
        fitmodel1$dPointObs(data1.point = data1.point,
                            model.point = model.point,
                            theta = theta)
      
      # Update state of the p particle
      state.particles[[p]] <- model.point
      
    }
    
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike1<- margLogLike1+ log(mean(weight.particles))
  }
  
  # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
  margLogLike2 <- 0
  
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  
  # Initialise time variable
  current.time <- 0
  
  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data2))){
    
    # Extract next data2 point (must be a vector)
    data2.point <- unlist(data2[i, ])
    next.time <- data2.point["time"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    
    ## Loop over particles: propagate and weight
    for(p in 1:n.particles){
      
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel2$simulate`
      traj <- fitmodel2$simulate(theta = theta,
                                 init.state = current.state.particle,
                                 times = c(current.time,next.time))
      
      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector
      model.point <- unlist(traj[2,fitmodel2$state.names])
      
      # Weight the particle with the likelihood of the observed 
      # data2 point using the function `fitmodel2$dPointObs`
      weight.particles[p] <-
        fitmodel2$dPointObs(data2.point = data2.point,
                            model.point = model.point,
                            theta = theta)
      
      # Update state of the p particle
      state.particles[[p]] <- model.point
      
    }
    
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike2<- margLogLike2 + log(mean(weight.particles))
  }
  
  # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
  margLogLike3 <- 0
  
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  
  # Initialise time variable
  current.time <- 0
  
  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data3))){
    
    # Extract next data3 point (must be a vector)
    data3.point <- unlist(data3[i, ])
    next.time <- data3.point["time"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    
    ## Loop over particles: propagate and weight
    for(p in 1:n.particles){
      
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel3$simulate`
      traj <- fitmodel3$simulate(theta = theta,
                                 init.state = current.state.particle,
                                 times = c(current.time,next.time))
      
      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector
      model.point <- unlist(traj[2,fitmodel3$state.names])
      
      # Weight the particle with the likelihood of the observed 
      # data3 point using the function `fitmodel3$dPointObs`
      weight.particles[p] <-
        fitmodel3$dPointObs(data3.point = data3.point,
                            model.point = model.point,
                            theta = theta)
      
      # Update state of the p particle
      state.particles[[p]] <- model.point
      
    }
    
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike3<- margLogLike3+ log(mean(weight.particles))
  }
  
  
  
  margLogLike <- margLogLike1 + margLogLike2 + margLogLike3
  
  ## Return marginal log-likelihood
  return(margLogLike)
  
}




my_mcmcMH <- function(target, init.theta, proposal.sd, n.iterations) {
  
  # evaluate the function "target" at "init.theta", and assign to
  # a variable called target.theta.current.
  target.theta.current <- target(init.theta)
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called theta.proposed.  
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta 
    # which will be useful when we will sample from a multivariate
    # target distribution
    
    theta.proposed <- rnorm(n = length(theta.current),mean = theta.current,sd = proposal.sd)
    
    # Note that 'rnorm' returns an unnamed vector, but the functions of
    # 'fitmodel' need a named parameter vector. We therefore set
    # the names of theta.proposed to be the same as the names of
    # theta.current
    names(theta.proposed) <- names(theta.current)
    
    # evaluate the function target at the proposed theta and
    # assign to a variable called target.theta.proposed
    target.theta.proposed <- target(theta.proposed)
    
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    log.acceptance <- target.theta.proposed - target.theta.current
    
    # draw random number number between 0 and 1 using "runif" and assign to
    # a variable called r.
    r <- runif(1)
    
    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    if (r < exp(log.acceptance)) {
      
      # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    
    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate 
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, theta.current, deparse.level=0)
    
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration)
    
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}


dTrajObs_POPDYN <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["N"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}


dLogPosPOPDYN <- function(fitmodel1, theta, init.state, data1) {
  log.prior <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood <- dTrajObs_POPDYN(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

dLogPosherpes <- function(fitmodel1,fitmodel2, theta, init.state, data1, data2) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_DE_hf(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_DE_hm(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
  
}

dTrajObs_POPDYN_sj <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["SJ"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}

dTrajObs_POPDYN_a <- function (fitmodel, theta, init.state, data, log = FALSE) {
  times <- rep(0:tail(data$time,n=1))
  traj <- fitmodel$simulate(theta, init.state, times)
  dens <- 0
  dens1 <- matrix(0,nrow(data),3)
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[data.point[2]+1, ])
    density <- fitmodel$dPointObs(data.point = data.point, 
                                  model.point = model.point, theta = theta, log = TRUE)
    dens1[i,] <-c(data.point[["obs"]],model.point[["A"]],density) 
    dens <- sum(dens1[,3])
  }
  return(ifelse(log, dens, exp(dens)))
}


dLogPosPOPDYN2 <- function(fitmodel1, fitmodel2, theta, init.state, data1, data2) {
  log.prior1 <- fitmodel1$dprior(theta, log = TRUE)
  log.likelihood1 <- dTrajObs_POPDYN_sj(fitmodel1, theta, init.state, data1, log = TRUE)
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- dTrajObs_POPDYN_a(fitmodel2, theta, init.state, data2, log = TRUE)
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
  
}

dLogPos_pmcmc <- function (fitmodel1,fitmodel2, theta, init.state, data1, data2, margLogLike = dTrajObs, 
          ...) 
{
  log.prior <- fitmodel1$dprior(theta = theta, log = TRUE)
  if (is.finite(log.prior)) {
    log.likelihood1 <- margLogLike(fitmodel = fitmodel1, theta = theta, 
                                  init.state = init.state, data = data1, ...)
  }
  else {
    log.likelihood <- -Inf
  }
  log.posterior1 <- log.prior1 + log.likelihood1
  
  log.likelihood2 <- margLogLike(fitmodel = fitmodel2, theta = theta, 
                                 init.state = init.state, data = data2, ...)
  
  log.posterior <- log.posterior1 + log.likelihood2
  
  return(log.posterior)
}

my_mcmcMH <- function(target, init.theta, proposal.sd, n.iterations) {
  
  # evaluate the function "target" at "init.theta", and assign to
  # a variable called target.theta.current.
  target.theta.current <- target(init.theta)
  
  # initialise variables to store the current value of theta, the
  # vector of samples, and the number of accepted runs
  theta.current <- init.theta
  samples <- theta.current
  accepted <- 0
  
  # run MCMC for n.iteration interations
  for (i.iteration in seq_len(n.iterations)) {
    
    # draw a new theta from the (Gaussian) proposal distribution
    # and assign to a variable called theta.proposed.  
    # See "?rnorm for more information
    # Note that this step is vectorized for any arbitratry theta 
    # which will be useful when we will sample from a multivariate
    # target distribution
    
    theta.proposed <- rnorm(n = length(theta.current),mean = theta.current,sd = proposal.sd)

    # Note that 'rnorm' returns an unnamed vector, but the functions of
    # 'fitmodel' need a named parameter vector. We therefore set
    # the names of theta.proposed to be the same as the names of
    # theta.current
    names(theta.proposed) <- names(theta.current)
    
    # evaluate the function target at the proposed theta and
    # assign to a variable called target.theta.proposed
    target.theta.proposed <- target(theta.proposed)
    
    # compute Metropolis-Hastings ratio (acceptance probability). Since
    # the multivariate Gaussian is symmetric, we don't need to consider
    # the proposal distribution here
    log.acceptance <- target.theta.proposed - target.theta.current
    
    # draw random number number between 0 and 1 using "runif" and assign to
    # a variable called r.
    r <- runif(1)
    
    # test acceptance by comparing the random number to the
    # Metropolis-Hastings ratio (acceptance probability) (using
    # "exp" because we calculated the logarithm of the
    # Metropolis-Hastings ratio before)
    if (r < exp(log.acceptance)) {
      
      # if accepted:
      # change the current value of theta to the proposed theta
      theta.current <- theta.proposed
      
      # updated the current value of the target
      target.theta.current <- target.theta.proposed
      
      # update number of accepted proposals
      accepted <- accepted + 1
    }
    
    # add the current theta to the vector of samples
    # Note that we use `rbind` in order to deal with multivariate 
    # target. So if `theta` is a vector then `samples` is a matrix.
    samples <- rbind(samples, theta.current, deparse.level=0)
    
    # print current state of chain and acceptance rate
    # use paste() to deal with the case where `theta` is a vector
    message("iteration: ", i.iteration, ", chain:", paste(theta.current, collapse=" "),
            ", acceptance rate:", accepted / i.iteration)
    
  }
  
  # return the trace of the chain (i.e., the vector of samples)
  return(samples)
}



my_particleFilter <- function(fitmodel, theta, init.state, data, n.particles) {
  
  ## Initialisation of the algorithm
  
  # Marginal log-likelihood is set to 0 and will be updated during the filtering steps
  margLogLike <- 0
  
  # Particle states can be stored in a list
  state.particles  <- rep(list(init.state), n.particles)
  
  # Weight: initially equal for all the particles 
  # particle weight can be stored in a vector
  weight.particles <- rep(1/n.particles, length = n.particles)
  
  # Initialise time variable
  current.time <- 0
  
  ## Loop over observation times: resample, propagate, weight
  for(i in seq_len(nrow(data))){
    
    # Extract next data point (must be a vector)
    data.point <- unlist(data[i, ])
    next.time <- data.point["time"]
    
    # Resample particles according to their weights. 
    # You can use the `sample` function of R
    # (normalisation of the weights is done in the function)
    index.resampled <- sample(x = n.particles,
                              size = n.particles,
                              replace = TRUE,
                              prob = weight.particles)
    state.particles <- state.particles[index.resampled]
    
    ## Loop over particles: propagate and weight
    for(p in 1:n.particles){
      
      # Extract current state of the particle 
      current.state.particle <- state.particles[[p]]
      
      # Propagate the particle from current observation time 
      # to the next one using the function `fitmodel$simulate`
      traj <- fitmodel$simulate(theta = theta,
                                init.state = current.state.particle,
                                times = c(current.time,next.time))
      
      # Extract state of the model at next observation time
      # Also make sure that model.point is a vector

      model.point <- unlist(traj[2,fitmodel$state.names])
      
      
      # Weight the particle with the likelihood of the observed 
      # data point using the function `fitmodel$dPointObs`
      weight.particles[p] <-
        fitmodel$dPointObs(data.point = data.point,
                           model.point = model.point,
                           theta = theta)
      
      # Update state of the p particle
      state.particles[[p]] <- model.point
      
    }
    
    # Increment time
    current.time <- next.time
    
    ## Increment the marginal log-likelihood
    # Add the log of the mean of the particles weights
    margLogLike <- margLogLike + log(mean(weight.particles))
  }
  
  ## Return marginal log-likelihood
  return(margLogLike)
  
}

postPredCheck_DE <- function(trace, n.samples, fitmodel, init.state, data) {
  
  # calculate maximum in obs column of data
  max.data <- max(data$obs)
  
  # draw n.samples random numbers between 1
  # and n.samples using the `samples` function 
  samples <- sample(seq_len(nrow(trace)), n.samples)
  
  # initialise vector of model maxima
  max.model <- c()
  
  # loop over samples
  for (i in samples) {
    
    # get i'th column from the trace, unlist
    # (to convert to a vector) and assign to parameter
    # vector theta
    theta <- unlist(trace[i, ])
    
    # use rObsTraj to generate
    # observation trajectory using theta
    times <- rep(0:tail(data$time,n=1))
    obs.traj <- rTrajObs(fitmodel, theta, init.state, times)
    
    # calculate maximum in model and add to max.model vector
    max.model <- c(max.model, max(obs.traj$obs))
  }
  
  # calculate quantiles of model maxima
  max.model.quant <- quantile(max.model, probs = c(0.025, 0.975))
  
  # calculate 2-sided p-value,
  # that is the proportion of elements of max.model which are
  # either greater or equal or less or equal (whichever is
  # less) and  multiply by 2 (because it is a 2-sided test)
  pvalue <- min(sum(max.model <= max.data),
                sum(max.model >= max.data)) / n.samples * 2
  
  # return two-sided p-value
  return(pvalue)
}

