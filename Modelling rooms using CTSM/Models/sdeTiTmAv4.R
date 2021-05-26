sdeTiTmAv4 <- function(data){
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  #room 1
  model$addSystem(dTi ~  1/Ci1*(1/Ria1*(Ta-T1) + 1/Rim1*(Tm1-T1)  + Ph1 + (a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5+a6*bs6+a7*bs7)*Gv)*dt + exp(p11)*dw1)
  model$addSystem(dTm1 ~  1/Cm1*(1/Rim1*(T1-Tm1))*dt + exp(p22)*dw2)
  #room 2
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-T2) + 1/Rim*(Tm2-T2)  + Ph1 + Aw*Gv)*dt + exp(p12)*dw3)
  model$addSystem(dTm2 ~  1/Cm*(1/Rim*(T2-Tm2))*dt + exp(p22)*dw4)
  #room 3
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-T3) + 1/Rim*(Tm3-T3)  + Ph2 + Aw*Gv)*dt + exp(p13)*dw5)
  model$addSystem(dTm3 ~  1/Cm*(1/Rim*(T3-Tm3))*dt + exp(p22)*dw6)
  #Room 4
  model$addSystem(dTi ~  1/Ci*(1/Ria*(Ta-T4) + 1/Rim*(Tm4-T4)  + Ph2 + (a1*bs1+a2*bs2+a3*bs3+a4*bs4+a5*bs5+a6*bs6+a7*bs7)*Gv)*dt + exp(p14)*dw7)
  model$addSystem(dTm4 ~  1/Cm*(1/Rim*(T4-Tm4))*dt + exp(p22)*dw8)
  
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph1,Ph2,bs1,bs2,bs3,bs4,bs5,bs6,bs7)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi1 ~ T1)
  model$addObs(yTi2 ~ T2)
  model$addObs(yTi3 ~ T3)
  model$addObs(yTi4 ~ T4)
  # Set the variance of the measurement error
  model$setVariance(yTi1 ~ exp(e11))
  model$setVariance(yTi2 ~ exp(e12))
  model$setVariance(yTi3 ~ exp(e13))
  model$setVariance(yTi4 ~ exp(e14))
  
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(T1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T4 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm4 = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e12 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e13 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(e14 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(a1 = c(init = 140, lb = -100, ub = 100))
  model$setParameter(a2 = c(init = 0, lb = 0.1, ub = 100))
  model$setParameter(a3 = c(init = 0, lb = 0.1, ub = 100))
  model$setParameter(a4 = c(init = 0, lb = 0.1, ub = 100))
  model$setParameter(a5 = c(init = 0, lb = 0.1, ub = 100))
  model$setParameter(a6 = c(init = 0, lb = 0.1, ub = 100))
  model$setParameter(a7 = c(init = 0, lb = 0.1, ub = 100))
  ##----------------------------------------------------------------    
  
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}