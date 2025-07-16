simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    S<-x[1]
    I<-x[2]
    E<-x[3]
    L<-x[4]
    Ps<-x[5]
    Pi<-x[6]
    Ldet <- x[7]

    
    with(as.list(params),{
      N <- S+I  # total host population
      P <- Ps + Pi # total parasite population
      b <- abs(0.537*sin(2*pi*(t/52-0.0587)))+0.537*sin(2*pi*(t/52-0.0587))
      K <- 62.47
      mu <- 0.05379
      alpha <- 1
      epsilon <- 0
      k <- 0.1
      tauh <- 42
      taue <- 717
      phih <- 1/8
      phie <- 1/9



      
      dS <- b*N*(K-N)/K + (gamma0 + gamma*Pi/I)*I - alpha*delta*S*E - mu*S - v*Ps
      dI <- alpha*delta*S*E - (mu + sigma+gamma0 + gamma*Pi/I)*I - v*Pi
      dE <- taue*I - phie*E
      dPs <- beta*L*S + (gamma0 + gamma*Pi/I)*Pi - alpha*delta*E*Ps - (mu + v + epsilon)*Ps - (Ps^2*v*(k+1))/(S*k)
      dPi <- beta*L*I + alpha*delta*E*Ps - (mu + v + epsilon + sigma + gamma0 + gamma*Pi/I)*Pi - (Pi^2*v*(k+1))/(I*k)
      dL <- tauh*P - phih*L
      dLdet <- tauh*P

      dx <- c(dS,dI,dE,dL,dPs,dPi,dLdet) 
    list(dx)
    })
  }
  
  # put incidence at 0 in init.state
  init.state["Ldet"] <- 0
  
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  
  traj$Ldet <- c(0,diff(traj$Ldet))
  
  return(traj)
  }

rPointObs <- function(model.point, theta){
  ##only model
  obs.point <- rnbinom(n=1,size=1,mu=model.point[["Ldet"]])
  return(c(obs=obs.point))
}

dPointObs <- function(data.point, model.point, theta, log = FALSE){
  ##model against data
  return(dnbinom(x=data.point[["obs"]],size=1,mu=model.point[["Ldet"]],log=log))
}

dprior <- function(theta, log = FALSE) {
  
  log.prior.delta <- dunif(theta[["delta"]], min = 1e-7, max = 1e-4, log = TRUE)
  log.prior.gamma0 <- dpois(round(1/theta[["gamma0"]]),lambda = 2.8, log = TRUE)
  log.prior.gamma <- dpois(round(1/theta[["gamma"]]),lambda = 10, log = TRUE)
  log.prior.beta <- dunif(theta[["beta"]], min = 1e-6, max = 1e-1, log = TRUE)
  log.prior.v <- dunif(theta[["v"]], min = 0, max = 1, log = TRUE)
  log.prior.sigma <- dunif(theta[["sigma"]], min = 0, max = 1, log = TRUE)
  
  
  log.sum = log.prior.delta + log.prior.gamma0+log.prior.beta + log.prior.v + log.prior.gamma + log.prior.sigma
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Coinfection model"
state.names <- c("S","I","E","L","Ps","Pi")
theta.names <- c("delta","gamma0","beta","v","gamma","sigma")

Coinf_det_w2 <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
