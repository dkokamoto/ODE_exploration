###########################################################################
###########################################################################
###                                                                     ###
### Example analysis of odes using exponential decay and seir models    ###
### Author: Examples McGillicutty                                       ###
### Last Update: 9/19/2022                                              ###
###                                                                     ###
###########################################################################
###########################################################################

# load required packages
library(RColorBrewer)
library(deSolve)

# load required functions 
source("R/example_odes.R")

### set parameters up here
scale = 365
maxtime = 3 * scale

# check out the ode function in deSolve package
?ode

### expontential decay example
# First, try the exponential decay function.  
# Lets set parameters on the annual scale but then scale to daily rates;
# this way we can visualize dynamics at small temporal scales

pars = c(z = 1 / scale)
init = c(x = 100) 

# generate ode solution at time = 0
exp_decay_ode(t = 0,
              state = init,
              parameters = pars)

# generate the solutions for the population sizes over time 
output <- ode(
  # supply initial values
  init, 
  # supply times you're interested in solving for
  times = seq(0, maxtime,by= 1), 
  # supply the function
  func = exp_decay_ode, 
  # supply the parameters
  pars)

# plot the numerical solution
plot(output[, -1], type = "l", ylab= "x", xlab = "time (days)", lwd = 4)

# calculate the analytical solution where N(t) = N(0) * exp(- z * t)
x_t = with(pars, init * exp(- pars * seq(0, maxtime,by= 1)))

# plot the analytical solution in red
lines(seq(0, maxtime,by= 1), x_t, lty= 3, col= "red", lwd= 2)

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 2 / scale)

# set some initial values (needs to be one for each state)
init = c(S= 100,
         E = 0, 
         I = 1,
         R = 0)

# generate the rate of change at time = 0
seir_ode(t= 0, 
         state = init, 
         pars
)

# generate the solutions for the population sizes over time 
output <- ode(
  # supply initial values
  init, 
  # supply times you're interested in solving for
  times = seq(0, maxtime,by= 1), 
  # supply the function
  func = seir_ode, 
  # supply the parameters
  pars)

# plot the output (the first column is the time column, so remove that one)
matplot(output[,- 1],type= "l", 
        col= brewer.pal(4,"RdYlBu"), 
        lty = 1,
        lwd = 2,
        ylab = "density",
        xlab = "time (days)")

# plot the legend so we can see whats happening
legend(1,
       max(output[,-1]), 
       names(init), 
       col= brewer.pal(4,"RdYlBu"), 
       xjust= -1, 
       bg = NA, 
       bty="n", 
       lty= 1)  

###############################################################################
### Logistic growth with threshold
###############################################################################
pars = c(z = 1 / scale,
         k = 2 / scale,
         r = 3 / scale)

init = c(x = 100) 

# generate the solutions for the population sizes over time 
output <- ode(
  # supply initial values
  init, 
  # supply times you're interested in solving for
  times = seq(0, maxtime,by= 1), 
  # supply the function
  func = log_growth_carry_cap_ode, 
  # supply the parameters
  pars)

library(dplyr)
library(ggplot2)
df <- as.data.frame(output) %>% mutate(individuals = x)

ggplot(df, aes(x = time, y = individuals))+
  geom_point()+
  labs(
    x = "Time (days)",
    y = "Number of Individuals",
    title = "Carrying Capacity"
  ) +
  gg_options()

###############################################################################
### Titration curves
###############################################################################
# parameters 
conc.base = 0.1
conc.acid = 0.1
pka1 = 5 
pka2 = 9
pkw = 14
vol.base = 50

# inflection points
inf_1 <- conc.base * vol.base/conc.acid
inf_2 <- 2*inf_1

# calculate theoretical titration
output <- make_titration(conc.base, 
                         conc.acid, 
                         pka1, 
                         pka2, 
                         pkw, 
                         vol.base)

# plot theoretical titration with inflection points
ggplot(output, aes(x = volume, y = ph))+
  geom_point()+
  geom_line(color = "blue")+
  geom_vline(xintercept = c(inf_1,inf_2), color = "red", linetype = 2)+
  labs(
    x = "Titrant Volume (mL)",
    y = "pH",
    title = "Theoretical Alkalinity Titration"
  )+
  gg_options()