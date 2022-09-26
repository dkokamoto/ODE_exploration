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
### Logistic growth with threshold 
###########################################################################
### parameters to be defined in `parameters`
# z = natural mortality

# r = growth rate

# k = constant of proportionality

### state variables to be defined in `state`
# x individuals left alive


log_growth_carry_cap_ode <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dx <- - z * x + (r * x)^2 - (k * x)^3
    list(c(dx))
  })
}

###########################################################################
### Theoretical titration curve for a diprotic weak base analyte using a 
### monoprotic strong acid as the titrant.
###########################################################################
make_titration <- function(conc.base, conc.acid, pka1, pka2, pkw, vol.base) {
  veq1 = conc.base * vol.base/conc.acid
  ka1 = 10^-pka1
  ka2 = 10^-pka2
  kw = 10^-pkw
  kw = 10^-pkw
  ph = seq(pkw, 1, -0.01)
  h = 10^-ph
  oh = kw/h
  delta = h - oh
  alpha1 = (ka1 * h)/(ka1 * ka2 + ka1 * h + h^2)
  alpha2 = h^2/(ka1 * ka2 + ka1 * h + h^2)
  volume = vol.base * 
    (conc.base * alpha1 + 2 * conc.base * alpha2 + delta)/
    (conc.acid - delta)
  df = data.frame(volume, ph)
  df <- df[df$volume > 0 & df$volume < 4 * veq1, ]
  rownames(df) = 1:nrow(df)
  return(df)
}

# Function to set up consistent ggplot aesthetics:
gg_options <-  function(){
  theme_bw(base_size = 16)+theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.background =  element_blank(),
    legend.background =  element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank())} 