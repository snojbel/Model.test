
rm(list=ls())  #removes stuff from environment.
# plot resource distribution            ----
par(mfrow=c(1,2))
z <- seq(0,20,0.01); omega<-7
plot( x=z, y=exp(- (z-9)^2 / (2 * omega)), lwd=4, col="gray", las=1, ylim=c(0,1), xlab="consumer trait, t", ylab=expression(paste("feeding ability, ",alpha)), type="l")
lines(x=z, y=exp(- (z-11)^2 / (2 * omega)), lwd=4, col="gray40")
text(x=2, y=0.9*1, "a)")

ct <- 8:12
plot(x=ct, y=rep.int(0.2, 5), las=1, ylab=expression(paste("frequency, ", p[r])), xlab=expression(paste("resource characteristic, ", c[i][,][r])), type="h", lwd=6, xlim=c(7,13), ylim=c(0,0.3))
text(x=8, y=0.30*0.9, "b)")

par(mfrow=c(1,1))

# define simulation function            ----
resourceCompetition <- function(resProp, resFreq, popSize, resGen=1, mutProb=0.001, mutVar=0.1, years=200, iniPmean=5, iniPvar=0.05, dispProb=0.5){

  # initialize population ......................................
  pop           <- matrix(NA, ncol=3, nrow=sum(popSize))      #NA is used here to preallocate memory for each entry while still not filling them with a value. 
  colnames(pop) <- c("patch", "pheno", "fec")
  
  pop[,1] <- c(rep.int(1,popSize[1]), rep.int(2, popSize[2]))
  pop[,2] <- rnorm(n=sum(popSize), mean=iniPmean, sd=0.05)
  
  stats         <- NULL
  phenotype <- matrix(NA, ncol=3, nrow=0)   # creates a matrix to store all phenotype values of all individuals in.
  colnames(phenotype) <- c("year", "patch", "phenotype")
  
  for(t in 1:years){
    
    # compute fecundity proxy ----
    # compute alpha - patch 1 ............................................................
    patch1    <- pop[pop[,1]==1,]   #This extracts all values in the pop matrix where the patch=1, i.e. they belong to patch 1. So all the rows with patch =2 are removed
    alphaSum1 <- NULL
    alpha1    <- NULL
    for(r in 1:ncol(resProp)){ # r<-1 # loop through each resource
      
      rp <- resProp[1,r]                                 # resource property 
      
      alpha     <- exp(-(patch1[,2]-rp)^2/ (2*resGen[1,1]) )  # compute all alphas
      alpha1    <- cbind(alpha1, alpha)                  # store all alphas
      alphaSum1 <- c(alphaSum1, sum(alpha) )             # store the sum over all alphas
      
    }
    
    # fec - patch 1 ............................................................
    for(i in 1:nrow(patch1)){
      
      fec <- 0
      for(r in 1:ncol(resProp)){ # r<-1; i<-1
        
        pi  <- resFreq[1,r]
        fec <- fec + pi * (alpha1[i,r]/alphaSum1[r])
        
      }
      patch1[i,3] <- fec
      
    }
    
    # alpha - patch 2 ............................................................
    patch2    <- pop[pop[,1]==2,]
    alphaSum2 <- NULL
    alpha2    <- NULL
    for(r in 1:ncol(resProp)){ # r<-1
      
      alpha     <- exp(-(patch2[,2]-resProp[2,r])^2/ (2*resGen[2,1]) )
      alpha2    <- cbind(alpha2, alpha)
      alphaSum2 <- c(alphaSum2, sum(alpha) )
      
    }
    
    # fec - patch 2 ............................................................
    for(i in 1:nrow(patch2)){
      
      fec <- 0
      for(r in 1:ncol(resProp)){ # r<-1; i<-1
        
        fec <- fec + resFreq[2,r]*alpha2[i,r]/alphaSum2[r]
        
      }
      patch2[i,3] <- fec
      
    }
    
    # create next generation  ----
    patch1new     <- patch1[sample(size=popSize[1], x=1:nrow(patch1),prob=patch1[,3], replace=T),]
    patch1new[,3] <- NA
    
    patch2new     <- patch2[sample(size=popSize[2], x=1:nrow(patch2),prob=patch2[,3], replace=T),]
    patch2new[,3] <- NA
    
    pop <- rbind(patch1new, patch2new)
    
    # mutate new generation   ----
    for(i in 1:nrow(pop)){
      
      if(runif(n=1, min=0, max=1)<=mutProb){
        
        pop[i,2] <- pop[i,2] + rnorm(n=1, mean=0, sd=mutVar)
        # print("mutation!")
        
      }
    }
    
    # dispersal               ----
    dispInds <- runif(n=nrow(pop))<dispProb
    dispInds <- which(dispInds==TRUE)
    if(length(dispInds)>=1){
      for(i in 1:length(dispInds)){
        if(pop[dispInds[i],1]==1){
          pop[dispInds[i],1] <- 2
        } else {
          pop[dispInds[i],1] <- 1
        }
      }
    }
    
    # extract stats           ----
    patch1    <- pop[pop[,1]==1,]
    patch2    <- pop[pop[,1]==2,]
    
    stats <- rbind( stats, c(t, mean(patch1[,2]), var(patch1[,2]), mean(patch2[,2]), var(patch2[,2]))) 
    #extract phenotypes of each individual each year
    phenotypes_patch1 <- cbind(rep(t, nrow(patch1)), rep(1, nrow(patch1)), patch1[,2])
    phenotypes_patch2 <- cbind(rep(t, nrow(patch2)), rep(2, nrow(patch2)), patch2[,2])
    
    phenotype <- rbind(phenotype, phenotypes_patch1, phenotypes_patch2)
  }
  
  # return output stats .............................................
  colnames(stats) <- c("year", "mean1", "var1", "mean2", "var2")
  colnames(phenotype) <- c("year", "patch", "phenotype")
  return(list(stats=stats, phenotype=phenotype))  #returns both the stats and the phenotype
  
}

# first test run                        ----
resFreqMatrix <- matrix(rep.int(0.2,5), nrow=2, ncol=5, byrow = TRUE); row.names(resFreqMatrix)<-c("patch1", "patch2")
resPropMatrix <- matrix(-2:2, nrow=2, ncol=5, byrow = TRUE)          ; row.names(resPropMatrix)<-c("patch1", "patch2")

output <- resourceCompetition(resProp=resPropMatrix, resFreq=resFreqMatrix, popSize=c(100, 100), resGen=matrix(c(0.2,0.2),ncol=1, nrow=2), mutProb=0.005, mutVar=0.5, years=250, iniPmean=1)

data <- output$stats

#For plotting "stats":

par(mfrow=c(2,1))
plot( x=data[,1], y=data[,2], lwd=2, col="skyblue", type="l", ylim=c(-2,2), las=1, xlab="year", ylab="avg. phenotype")
lines(x=data[,1], y=data[,4], lwd=2, col="coral"); abline(h=0, col="gray", lty=1)
legend("top", legend=c("patch 1", "patch 2"), lwd=2, lty=1, col=c("skyblue", "coral"), box.col="transparent", horiz=TRUE)

plot( x=data[,1], y=data[,3], lwd=2, col="skyblue", type="l", ylim=c(0,4), las=1, xlab="year", ylab="phenotypic variance")
lines(x=data[,1], y=data[,5], lwd=2, col="coral")
par(mfrow=c(1,1))

#Plotting "phenotypes":

phenotype_data <- output$phenotype 


# Create a scatter plot of individual phenotypes
plot(phenotype_data[,1], phenotype_data[,3], pch=19, col=rgb(0.5,0.2,0.5, alpha =0.3), xlab="Year", ylab="Phenotype")

# Add a legend to distinguish between patches
legend("topright", legend=c("Patch 1", "Patch 2"), col=c(1, 2), pch=19)

# Add a title to the plot
title("Individual Phenotypes Over Years for Patch 1 and Patch 2")



# explore parameters with replicates    ----
propVal    <- c(0.25, 0.4, 0.50, 0.55, 0.6, 0.7, 1.00)
propVar <- NULL

for(i in 1:length(propVal)){ # i<-1
  propVar    <- c(propVar, sum(((-2:2)*propVal[i]-0)^2)/5 )
}

dataP_mean <- NULL
dataP_sd   <- NULL
for(p in 1:length(propVal)){
  dataR_mean <- NULL
  dataR_sd   <- NULL
  for(r in 1:2){ 
    data       <- resourceCompetition(resProp=resPropMatrix*propVal[p], 
                                      resFreq=resFreqMatrix, 
                                      resGen=matrix(c(0.5,0.5),ncol=1, nrow=2), 
                                      mutProb=0.005, mutVar=0.5, 
                                      popSize=c(1000, 1000), 
                                      years=250, iniPmean=1)[250, ,drop=F]
    dataR_mean <- c( dataR_mean, mean(c(data[,2], data[,4])) )
    dataR_sd   <- c( dataR_sd  , mean(c(data[,3], data[,5])) )
  }
  dataP_mean <- rbind(dataP_mean, c(propVar[p], mean(dataR_mean), max(dataR_mean), min(dataR_mean)))
  dataP_sd   <- rbind(dataP_sd  , c(propVar[p], mean(dataR_sd)  , max(dataR_sd)  , min(dataR_sd)))
  print(p)
  
}
colnames(dataP_mean) <- c("resource var", "mean z", "max z"     , "min z")
colnames(dataP_sd)   <- c("resource var", "mean z", "max var(z)", "min var(z)")

# plot  .......................
par(mfrow=c(2,1))
plot(x=dataP_mean[,1], y=dataP_mean[,2], lwd=2, type="o", las=1, xlab="resource variation", ylab="mean phenotype", ylim=c(-1,1), xlim=c(0,max(dataP_sd[,1])))
lines(x=dataP_mean[,1], y=dataP_mean[,3], lwd=2, col="gray")
lines(x=dataP_mean[,1], y=dataP_mean[,4], lwd=2, col="gray")
abline(h=0, col="gray", lty=3)

plot(x=dataP_sd[,1], y=dataP_sd[,2], lwd=2, type="o", las=1, xlab="resource variation", ylab="phenotypic variance", ylim=c(0,2), xlim=c(0,max(dataP_sd[,1])))
lines(x=dataP_sd[,1], y=dataP_sd[,3], lwd=2, col="gray")
lines(x=dataP_sd[,1], y=dataP_sd[,4], lwd=2, col="gray")
abline(v=0.5, col="dodgerblue", lwd=2, lty=3)
legend("top", legend=c("resource generalism"), lty=3, lwd=2,col="dodgerblue", box.col="transparent")

