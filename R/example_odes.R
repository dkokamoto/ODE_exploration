###########################################################################
###########################################################################
###                                                                     ###
### Sets of ODE models for analysis                                     ###
### Author: Examples McGillicutty                                       ###
### Last Update: 9/19/2022                                              ###
###                                                                     ###
###########################################################################
###########################################################################

### ode functions take three arguments:
# '@param t vector of time  points
# '@param state vector of state variables
# '@param parameters vector of parameter values 


###########################################################################
### simple exponential decay 
###########################################################################

# This model is simply exponential decay over time.
# Note: The analytical solution is x(t) = x(0) * exp(- z * t)

### parameters to be defined in `parameters`
# z = natural mortality

### state variables to be defined in `state`
# x individuals left alive

exp_decay_ode <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dx <- - z * x 
    list(c(dx))
  })
}

###########################################################################
### SEIR model with fishing
###########################################################################

# This model characterizes changes in susceptible, exposed, infected, 
# and resistant/recovered individuals over time as a function of baseline
# mortality rate (note that one can convert this into annual surival by 
# exp(- z * t)) where t = the time in a year relative to the scale of z 


### parameters to be defined in `parameters`
# z = instantaneous natural mortality (per individual, per unit time)
# f = instantaneous fisheries mortality (per individual, per unit time)
# beta = instantaneous rate of infection (per infected, per susceptible, per unit time)
# k = instantaneous rate of transition from exposed to infected
# r = instantaneous rate of recruitment into susceptible
# delta = instantaneous death rate of infected over and above natural/fisheries 
# gamma = instantaneous rate of transition into recovered from infected 

### state variables to be defined in `state`
# S = susceptible
# E = exposed
# I = infected/infectious
# R = recovered/resistant

seir_ode <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # susceptible
    dS <- - (z + f + beta * I) * S  + r 
    # exposed but not yet infectious
    dE <- - (z + f + k) * E + beta * I * S 
    # infected
    dI <- k * E - (z + f + delta + gamma) * I
    # recovered/resistant
    dR <- gamma * I - (z + f) * R
    # return the rate of change for the compartments
    list(c(dS, dE, dI, dR))
  })
}

###########################################################################
### SEIR model without fishing
###########################################################################

# This model characterizes changes in susceptible, exposed, infected, 
# and resistant/recovered individuals over time as a function of baseline
# mortality rate (note that one can convert this into annual surival by 
# exp(- z * t)) where t = the time in a year relative to the scale of z 


### parameters to be defined in `parameters`
# z = instantaneous natural mortality (per individual, per unit time)
# beta = instantaneous rate of infection (per infected, per susceptible, per unit time)
# k = instantaneous rate of transition from exposed to infected
# r = instantaneous rate of recruitment into susceptible
# delta = instantaneous death rate of infected over and above natural/fisheries 
# gamma = instantaneous rate of transition into recovered from infected 

### state variables to be defined in `state`
# S = susceptible
# E = exposed
# I = infected/infectious
# R = recovered/resistant

seir_nofish_ode <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # susceptible
    dS <- - (z + beta * I) * S  + r 
    # exposed but not yet infectious
    dE <- - (z + k) * E + beta * I * S 
    # infected
    dI <- k * E - (z + delta + gamma) * I
    # recovered/resistant
    dR <- gamma * I - (z) * R
    # return the rate of change for the compartments
    list(c(dS, dE, dI, dR))
  })
}

