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
maxtime = 15 * scale

# check out the ode function in deSolve package
?ode

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 2 / scale,
          zeta=0,
          theta=0,
          f_sick=0)

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

## WITH SOME FATAL INFECTIONS

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 1 / scale,
          zeta= 1 / scale,
          theta = 0, 
          f_sick=0)

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
        xlab = "time (days)",
        main = "Disease mort")

# plot the legend so we can see whats happening
legend(1,
       max(output[,-1]), 
       names(init), 
       col= brewer.pal(4,"RdYlBu"), 
       xjust= -1, 
       bg = NA, 
       bty="n", 
       lty= 1)

## WITH SOME FATAL INFECTIONS AND SOME IMMUNITY LOSS

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 1 / scale,
          zeta= 1 / scale,
          theta = 0.5 / scale, 
          f_sick=0)

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
        xlab = "time (days)", 
        main="Disease mort + immunity loss")

# plot the legend so we can see whats happening
legend(1,
       max(output[,-1]), 
       names(init), 
       col= brewer.pal(4,"RdYlBu"), 
       xjust= -1, 
       bg = NA, 
       bty="n", 
       lty= 1)

## WITH SOME FATAL INFECTIONS AND SOME IMMUNITY LOSS, AND HARVEST

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0.1 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 1 / scale,
          zeta= 1 / scale,
          theta = 0.5 / scale, 
          f_sick=0.1)

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
        xlab = "time (days)",
        main="All oysters harvest f=0.1 + disease\nmort + immunity loss")

# plot the legend so we can see whats happening
legend(1,
       max(output[,-1]), 
       names(init), 
       col= brewer.pal(4,"RdYlBu"), 
       xjust= -1, 
       bg = NA, 
       bty="n", 
       lty= 1)

## WITH SOME FATAL INFECTIONS AND SOME IMMUNITY LOSS, AND SELECTIVE HARVEST

# set parameters
pars =  c(z = 0.1 / scale,  
          f = 0.2 / scale, 
          beta = 2 / scale, 
          k = 2 / scale,
          r = 10 / scale, 
          delta = 0.1 / scale,
          gamma = 1 / scale,
          zeta= 1 / scale,
          theta = 0.5 / scale, 
          f_sick=0)

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
        xlab = "time (days)",
        main="Healthy oysters harvest f=0.2 + disease\nmort + immunity loss")

# plot the legend so we can see whats happening
legend(1,
       max(output[,-1]), 
       names(init), 
       col= brewer.pal(4,"RdYlBu"), 
       xjust= -1, 
       bg = NA, 
       bty="n", 
       lty= 1)
