sdeTiTm6 <- function(data, yTi,Ph){
  data$yTi <- yTi
  data$Ph <- Ph
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  
 
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-Ti) + 1/Rim*(Tm-Ti)+1/Rie*(Te-Ti) 
                               + (1/Rih)*(Th-Ti) + 
                                 p*(a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5)*Gv)*dt + exp(p11)*dw1)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti-Tm))*dt + exp(p22)*dw2)
  model$addSystem(dTh ~ 1/Ch*(1/Rih*(Ti-Th) + Ph)*dt + exp(p33)*dw3 )
  model$addSystem(dTe ~ 1/Ce*(1/Rie*(Ti-Te) + 1/Rea*(Ta-Te) + (1-p)*(a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5)*Gv)*dt + exp(p44)*dw4 )
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph, bs1, bs2,bs3,bs4,bs5)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi ~ Ti)
  # Set the variance of the measurement error
  model$setVariance(yTi ~ exp(e11))
  
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti = c(init = 20, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 20, lb = 0, ub = 40))
  model$setParameter(Th = c(init=15  ,lb=-15, ub=80) )
  model$setParameter(Te = c(init = 10, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 250, lb = 1E-5, ub = 1E5))
  model$setParameter(Ch = c(init=20   ,lb=1E-5  ,ub=1E5))
  model$setParameter(Ce = c(init=10  ,lb=1E-5  ,ub=1E5))
  model$setParameter(p = c(init=0.5  ,lb=1E-5  ,ub=0.9999))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rih = c(init = 1   ,lb=1E-3  ,ub=1E3))
  model$setParameter(Rie = c(init = 1   ,lb=1E-3  ,ub=1E3))
  model$setParameter(Rea = c(init = 1   ,lb=1E-3  ,ub=1E3))
  
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 20))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 20))
  model$setParameter(p33 = c(init=  -1  ,lb=-30   ,ub= 20) )
  model$setParameter(p44 = c(init=  -1  ,lb=-40   ,ub= 20) )
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  
  model$setParameter(a1 = c(init = 85, lb = 1E-2, ub = 140))
  model$setParameter(a2 = c(init = 6, lb = 1E-2, ub = 30))
  model$setParameter(a3 = c(init = 6, lb = 1E-2, ub = 30))
  model$setParameter(a4 = c(init = 3, lb = 1E-2, ub = 30))
  model$setParameter(a5 = c(init = 3, lb = 1E-2, ub = 30))
  ##----------------------------------------------------------------    
  
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}