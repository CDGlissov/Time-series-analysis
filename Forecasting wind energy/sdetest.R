windpow2 <- function(data){
  data$dWs1 = c(0,diff(data$Ws1))
  # Generate a new object of class ctsm
  model = ctsm()
  model$options$odeeps <- 1E-4
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dX ~  ( ((1-exp(-X)) * (rho*dWs1+R)+theta_x*(Ws1*mu-X))/(sigma_x*sqrt(X)) - sigma_x/(4*sqrt(X)) )*dt + dw1)
  model$addSystem(dR ~  -theta_r*R*dt + sigma_r*dw2)
  # Set the names of the inputs
  model$addInput(Ws1,dWs1)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(p ~ (0.5-0.5*(sinh(gamma2*(X-gamma3))/cosh(gamma2*(X-gamma3))))*z3/(1+exp(-z1*(X-z2))) )
  # Set the variance of the measurement error
  model$setVariance(p ~ exp(e11))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(X = c(init = 5, lb = 1E-4, ub = 25))
  model$setParameter(R = c(init = 1, lb = 1E-4, ub = 20))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  #model$setParameter(mu = c(init = 5, lb = 1E-4, ub = 10))
  model$setParameter(mu = c(init = 0.7120, lb = 1E-4, ub = 4))
  #model$setParameter(rho = c(init = 0.110921, lb = 1E-4, ub = 4))
  model$setParameter(rho = c(init = 0.3087, lb = 1E-4, ub = 1))
  model$setParameter(sigma_x = c(init = 0.0934538, lb = 1E-4, ub = 1))
  model$setParameter(sigma_r = c(init = 0.9427243, lb = 1E-4, ub = 1))
  
  
  model$setParameter(theta_x = c(init = 0.5489318, lb = 1E-4, ub = 1))
  model$setParameter(theta_r = c(init = 1.358, lb = 1E-4, ub = 5))
  
  #model$setParameter(theta_r = c(init = 0.45, lb = 1E-4, ub = 1))
  
  #model$setParameter(e11 = c(init = 3.9964, lb = -1E+2, ub = 1E+1))
  
  model$setParameter(e11 = c(init = 4, lb = -1E+2, ub = 10))
  #model$setParameter(gamma2 = c(init = 1, lb = 1E-4, ub = 2))
  model$setParameter(gamma2 = c(init = 0.614, lb = 1E-4, ub = 1))
  model$setParameter(z1 = c(init = 0.612, lb = 1E-4, ub = 5))
  #model$setParameter(z1 = c(init = 0.9, lb = 1E-4, ub = 1))
  model$setParameter(z2 = c(init = 8, lb = 1E-4, ub = 10))
  
  
  #iteration 2
  model$setParameter(z3 = c(init = 0.5, lb = 1E-4, ub = 1))
  model$setParameter(gamma3 = c(init = 19, lb = 1E-4, ub = 40))
  ##----------------------------------------------------------------    
  
  # Run the parameter optimization
  fit = model$estimate(data=data)
  return(fit)
}

