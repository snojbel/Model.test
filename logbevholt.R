

rm(list=ls())
set.seed(321)


# THE DETERMINISTIC MODEL ................................ ----
# pop growth        ----
# function definition:
# this function computes the deterministic population growth IN ABSENCE OF STOCHASTICITY (click on small triangle to unfold the function definition)
logisticGrowth_BevHolt_math <- function(N0=2, b=0.001, l_max=1.4, years=100){ 
  
  Nadlt    <- N0                                 # set initial adult number (with input parameter "N0")
  dataMath <- NULL                               # set object where you want to save the summary statistics
  dataMath <- rbind(dataMath, c(0, N0))          # save summary stats for initial condition (1st column: year=0, 2nd column, N=N0)
  
  for(t in 1:years){                             # loop through the years (where the number of years is determined by input parameter "years")
    
    indFec   <- l_max / (1 + b * Nadlt)          # compute avg. individual fecundity following Beverton-Holt (with input parameters "b" and "l_max")
    Noff     <- Nadlt * indFec                   # compute total number of offspring individuals
    Nadlt    <- Noff                             # all adults die & offspring become adults (offspring overwrite adults)
    dataMath <- rbind(dataMath, c(t, Nadlt))     # save stats
    
  }
  colnames(dataMath) <- c("year", "N")
  
  return(dataMath)
  
}

# function usage:
dataMath <- logisticGrowth_BevHolt_math(N0=2, b=0.001, l_max=1.4, years=50)
plot( x=dataMath[,1], y=dataMath[,2], col="skyblue", lty=1, lwd=3, las=1, type="p", pch=16, xlab="time", ylab="population size (N)", ylim=c(0,500)) # 

# fecundity         ----
# in the following lines, the fecundity is plotted (y-axis) against population size (x-axis)
b     <- 0.001   # competition coefficient
l_max <- 1.4     # max. fecundity at very low densities
N     <- 1:1000  # possible population sizes
fec   <- l_max / (1 + b * N)
plot(x=N, y=fec, type="p", pch=16, lwd=3, col="coral", las=1, ylim=c(0,1.5), xlab="population size, N", ylab="individual fecundity, f(N)")
abline(h=1, col="gray", lty=2)
rm(b, l_max, N, fec) # removing the parameters from the workspace again.

# carrying capacity ----
# this function computes the expected carrying capacity from b and l_max
logisticGrowth_BevHolt_K <- function(b=0.001, l_max=1.4){ 
  # this function computes the expected carrying capacity IN ABSENCE OF STOCHASTICITY
  
  carryingCapacity <- (l_max - 1) / b
  return( carryingCapacity )
  
}
abline(v=logisticGrowth_BevHolt_K(b=0.001, l_max=1.4)) # adding a line to a plot using the function abline, while the parameter v stands for "vertical"

#                   ----
# INDIVIDUAL-BASED SIMULATIONS ........................... ----
# stochastic fec    ----
# we will sample the individual fecundity from a Poisson distribution (check wikipedia for more information)
# The Poisson is a discrete distribution (result in whole numbers / integers). This makes sense for fecundity.

hist(rpois(n=100000, lambda=2)+0.001, main="Poisson distribution", xlab="individual fecundity, f(N)", las=1, freq = F) 

# change lambda to see the probability distribution with other mean fecundities.

# define function   ----
# this function simulates one replicate of logistic growth with stochastic fecundity
logisticGrowth_BevHolt_sim <- function(N0=2, b=0.001, l_max=1.4, years=50){
  
  # example input ........................................................
  # inside a function, you can only use the input parameters (the parameters inside of function() ...). Every other parameter in the workspace cannot be used.
  # When writing your function, it is often helpful to define a set of example input parameters to see how the function works. Uncomment the following line (remove the hashtag) and run it.
  # N0<-2; b<-0.001; l_max<-1.4; years<-50 # example input parameters for coding purposes
  # However, comment the previous line again when you are done with coding the function !!! Otherwise, your input parameters will always be overwritten.
  
  # 0) set initial population matrix .....................................
  ADLT           <- matrix(1, ncol=1, nrow=N0)     # Create a population called "ADLT" that consists of N0 individuals (the initial population size) with one state variable "patch"
  colnames(ADLT) <- "patch"                        # Here, each individual is located in patch=1
  
  sumStats      <- NULL                            # Create the object "sumStats" that stores the summary statistics during your simulations. Let it be empty in the beginning.
  sumStats      <- rbind(sumStats, c(0,N0))        # Then, fill this object with the initial conditions (time is zero, and population size is N0)
  
  # loop through the annual events ......................................
  for(t in 1:years){ # t<-2
    
    # 1) reproduction (adults produce offspring)
    N       <- nrow(ADLT)                           # At the beginning of each year, first compute the population size "N"
    Fec     <- l_max / (1+b*N)                      # Then compute the average fecundity of each adult based on this population size N (and the input parameters b and l_max)
    indFec  <- rpois(n=N, lambda=Fec)               # STOCHASTIC INDIVIDUAL FECUNDITY: now pick for each of the N individuals the fecundity from a poisson distribution with mean=Fec
    OFF     <- matrix(1, ncol=1, nrow= sum(indFec)) # Create the offspring population (where the row number is identical to the sum of the adult fecundities).
    colnames(OFF) <- "patch"                        # Again, all offspring are located in patch=1.
    
    # 2) aging (adults die, offspring become adults)
    ADLT     <- OFF                                 # overwrite the matrix ADLT with OFF (that means that all adults die and the entire population is now made up of offspring that just became adults)
    rm(OFF)                                         # remove the matrix OFF from workspace (not obligatory, but good coding practice)
    
    # 3) store summary statistics
    sumStats <- rbind(sumStats, c(t,nrow(ADLT)))    # extract summery statistics and store them in sumStats
    
  }
  colnames(sumStats) <- c("year", "size")
  
  return(sumStats)
  
  # ! don't forget to put back the # before the example input parameter
  
}

# first test run    ----
# run a first test simulation and have a look at the output
simOutput <- logisticGrowth_BevHolt_sim(N0=10, b=0.001, l_max=1.4, years=50) 

plot(x=simOutput[,column=1], y=simOutput[,column=2], lwd=2, type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,500))

# debugging         ----
# 1) personal code checking (does the code make sense?)
# 2) error message in R ?
# 3) do some plots and play a bit with the parameters (do the results make sense?)
# 4) compare to mathematical expectation
# a) carrying capacity
simOutput <- logisticGrowth_BevHolt_sim(N0=10, b=0.001, l_max=1.4, years=100)
K         <- logisticGrowth_BevHolt_K(         b=0.001, l_max=1.4)
plot(  x=simOutput[,column=1], y=simOutput[,column=2], lwd=2, type="l", xlab="time", ylab="population size (N)", las=1)
abline(h=K, lty=3, col="gray", lwd=2) # this adds a line at the carrying capacity predicted by the mathematical model
rm(K)

# b) growth tragectory
dataMath <- logisticGrowth_BevHolt_math(N0=10, b=0.001, l_max=1.4, years=100)
lines(x=dataMath[,1], y=dataMath[,2], col="skyblue", lty=2, lwd=2)

# 5) extreme scenario tests
# a) what happens when we start with zero individuals ?
simOutput <- logisticGrowth_BevHolt_sim(N0=0, b=0.001, l_max=1.4, years=100)
plot(x=simOutput[,1], y=simOutput[,2], lwd=2, type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,500))
abline(h=logisticGrowth_BevHolt_K(b=0.001, l_max=1.4), lty=3, col="gray", lwd=2) # this is the carrying capacity predicted by the mathematical model

# b) what happens when we start with with population size larger than the carrying capacity ?
simOutput <- logisticGrowth_BevHolt_sim(N0=600, b=0.001, l_max=1.4, years=100)
plot(x=simOutput[,1], y=simOutput[,2], lwd=2, type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,650))
abline(h=logisticGrowth_BevHolt_K(b=0.001, l_max=1.4), lty=3, col="gray", lwd=2) # this is the carrying capacity predicted by the mathematical model

# c) what happens at extreme fecundities ?
simOutput <- logisticGrowth_BevHolt_sim(N0=10, b=0.001, l_max=40, years=100)
plot(x=simOutput[,1], y=simOutput[,2], lwd=2, type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,50000))
abline(h=logisticGrowth_BevHolt_K(b=0.001, l_max=40), lty=3, col="gray", lwd=2) # this is the carrying capacity predicted by the mathematical model


#                   ----
# SIMULATION DESIGN ...................................... ----
# replicates        ----
set.seed(134)   # this allows to replicate results, despite stochasticity. Everybody running this script then should get the same results/plot
simOutput <- NULL
for(r in 1:10){ # running 10 replicates (independent simulation runs)
  
  simOutput <- cbind(simOutput, logisticGrowth_BevHolt_sim(N0=3, b=0.001, l_max=1.4, years=50)[,2])  # store the population sizes of each replicate in a single matrix
  
}
# plot the replicates beside each other:
plot( x=0:50, y=simOutput[,1], lwd=1, col="gray", type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,500))
lines(x=0:50, y=simOutput[,2], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,3], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,4], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,5], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,6], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,7], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,8], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,9], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,10], lwd=1, col="gray")

dataMath <- logisticGrowth_BevHolt_math(N0=3, b=0.001, l_max=1.4, years=50)
lines(x=dataMath[,1], y=dataMath[,2], col="black", lty=1, lwd=2)

# init. conditions  ----
# when using a slightly different seed, we "suddenly" observe population extinction. 
set.seed(111)    
simOutput <- NULL
for(r in 1:10){ # running 10 replicates (independent simulation runs)
  
  simOutput <- cbind(simOutput, logisticGrowth_BevHolt_sim(N0=3, b=0.001, l_max=1.4, years=50)[,2])  # store the population sizes of each replicate in a single matrix
  
}
# plot the replicates beside each other:
plot( x=0:50, y=simOutput[,1], lwd=1, col="gray", type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,500))
lines(x=0:50, y=simOutput[,2], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,3], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,4], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,5], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,6], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,7], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,8], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,9], lwd=1, col="gray")
lines(x=0:50, y=simOutput[,10], lwd=1, col="gray")

dataMath <- logisticGrowth_BevHolt_math(N0=3, b=0.001, l_max=1.4, years=50)
lines(x=dataMath[,1], y=dataMath[,2], col="black", lty=1, lwd=2)
# What is happening here?

# Instead, when we would start our simulations with much larger N0, we would never observe population extinction. Run the following lines over and over again to test this.
simOutput <- logisticGrowth_BevHolt_sim( N0=100, b=0.001, l_max=1.4, years=100)
dataMath  <- logisticGrowth_BevHolt_math(N0=100, b=0.001, l_max=1.4, years=100)
plot( x=0:100, y=simOutput[,2], lwd=1, col="gray", type="l", xlab="time", ylab="population size (N)", las=1, ylim=c(0,650))
lines(x=dataMath[,1], y=dataMath[,2], col="black", lty=1, lwd=2)
# This result illustrates that the initial simulation conditions matter. Chose them wisely.

# burn-in           ----
# When we run simulations, the simulated population often changes quite strongly in the beginning (as in the previous examples).
# These dynamics heavily depend on the initial conditions.
# However, in many occasions we want to come to conclusions that are independent of the initial conditions. 
# For instance, sometimes we want to know the equilibrium state when all processes balanced out.
# And even when we are interested in temporal dynamics, we need to initiate these dynamics from an equilibrium state. 
# Otherwise we cannot be sure whether the dynamics are the result of the initial condition or date back to some perturbation (for instance environmental changes).
# Therefore, you need to run the simulations for a "burn-in" period, and only after this "burn-in" extract summary statistics or initiate temporal changes.
# Here is an example
simOutput <- NULL
for(r in 1:200){  # simulation of 200 replicates
  # run simulations for 200 years of "burn-in" and store all N values in a matrix. In "simOutput", each column holds the N values of a single replicate. Each row holds data for a specific year
  simOutput <- cbind(simOutput, logisticGrowth_BevHolt_sim(N0=100, b=0.001, l_max=1.4, years=200)[,2]) 
}
Nt200 <- simOutput[200,]            # here we only use the data of the last year (row=200)
mean(Nt200); min(Nt200); max(Nt200) # and compute the avg. over replicates afte the burn-in
hist(Nt200, xlim=c(200,600))        # Recall that the math. model predicts a carrying capacity of 400
# Also note that this result changes quite dramatically when starting the simulations with a single individual (N0=1)
# And we get very different result when teh burn-in period is to short (e.g. after only 10 years). 

# parameter space   ----
lambda_val <- c(seq(1.5,5,0.5))                                        # here, we specify which l_max values we want to explore
summary    <- NULL                                                     # define object to collect data for each l_max value
for(a in 1:length(lambda_val)){                                        # loop through the parameter space (all l_max values)
  
  simOutput <- NULL                                                    # set summary data for specific l_max to NULL
  for(r in 1:100){                                                     # simulate 200 replicates for each l_max
    # Run simulations for 200 years of "burn-in" and add data to sinOutput.
    # Note that the "l_max" value is updated for every loop cycle of the outer loop ("a")! That means in the first loop cycle through the parameter space a=1 and l_max will be lambda_val[1], the first entry of lambda_val. In the second cycle, a=2 ...
    simOutput <- cbind(simOutput, logisticGrowth_BevHolt_sim(N0=100, b=0.001, l_max=lambda_val[a], years=200)[,2]) 
  }
  Nt200   <- simOutput[200,]                                             # here we only use the data of the last year (= 200 years of burn-in)
  summary <- rbind(summary, c(lambda_val[a], mean(Nt200), min(Nt200), max(Nt200)))  # Store the average over 100 replicates after 200 years of burn-in in summary
  
}
colnames(summary) <- c("l_max", "avg. N", "min N", "max N")

plot(x=summary[,1], y=summary[,2], pch=16, las=1, xlab="l_max", ylab="mean N", ylim=c(0,5000))
points(x=lambda_val, y=logisticGrowth_BevHolt_K(b=0.001, l_max=lambda_val), col="skyblue", cex=2)
segments(x0=summary[,1], y0=summary[,3], x1=summary[,1], y1=summary[,4], col="black", lwd=1)

#                                                          ----
#                   ----
#                   ----


