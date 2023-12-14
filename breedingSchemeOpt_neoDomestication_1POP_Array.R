#if (!require("AlphaSimR")) install.packages("AlphaSimR")
library(AlphaSimR)
#if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
#if (!require("sommer")) install.packages("sommer")
library(sommer)
#if (!require("optiSel")) install.packages("optiSel")
library(optiSel)
#if (!require("rrBLUP")) install.packages("rrBLUP")
library(rrBLUP)

# Breeding Scheme Design and Optimization for Neo-Domestication
### The goal of this code is to perform stochastic simulation of genetic gain in populations that vary in
### their domestication status (wild, orphan, landrace). The design will follow 2 parts: (1) Population
### formation (burn-in) phase; and (2) Breeding scheme development phase. 
### _______________________________________________________________________________________________________
### Phase 1 is to form populations typically expected within the above listed types by altering effective 
### population size (Ne), setting alternative germplasm collection counts (Nfounder), and setting different
### number of breeding parents (Nparents). Following this, 2 traits will be added to the founding population: 
### (1) simple oligogenic trait to represent disease tolerance, seed size and early vigor, or other; 
### and (2) complex polygenic trait to represent yield or adaptive plasticity. The covariation of traits 
### will vary from being negatively correlated (-0.5), to uncorrelated (0), to positively correlated (0.5).
### Following this, different forms and degrees of random or directional (phenotypic) selection will be 
### applied across 40 cycles: random selection for the wild population, low intensity of selection for the 
### orphan population, and low-medium intensity for the landrace population. 
### _______________________________________________________________________________________________________
### Phase 2 is to develop breeding schemes to identify parameters of interest that increase the speed to
### achieve high gain with low variation for trait 1 (simple) and high gain with high variation for trait 2
### (complex). This will be done by setting different number of parents (Nparents2), varying number of 
### populations formed (Npops) as well as the timing of their recombination (Popcomb), the number of 
### environments (Env) and replications within environments (Reps), and selecting using a weighted selection
### index with different weights applied to each trait, performing selection using different criteria 
### (Selcrit) and intensity (Selintens).
### _______________________________________________________________________________________________________
### The outputs of this stochastic simulation through the use of for-loops will include all parameters set,
### genetic gain realized, and genetic variance lost. Also, the cost of each scheme will be estimated. These
### outputs will then be used in mixed-model analysis to estimate parameters which maximize the return
### on investment (ROI). 

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assuming your parameters are space-separated values on a single line
# Split the arguments into separate objects
parameters <- unlist(strsplit(args, " "))

# Now, you can access and use the parameters as objects in your R script
Ne <- as.numeric(parameters[1])
nPa <- as.numeric(parameters[2])
nPr <- as.numeric(parameters[3])
h1 <- as.numeric(parameters[4])
h2 <- as.numeric(parameters[5])
sel <- parameters[6]
siC <- as.numeric(parameters[7])
h <- c(h1,h2)

# Phase 1: Population Formation (Burn-in)
### Set phase 1 parameters
z1QTL <- c(0,1,0,1,3,0,0,1,2,0) # simple oligogenic trait (trait 1: 8 QTL)
z2QTL <- c(40,10,40,100,30,80,20,10,20,50) # complex polygenic trait (trait 2: 400 QTL)
z1mean <- 0; z1var <- 1
z2mean <- 0; z2var <- 1
popType <- c("wild","orphan","landrace")

### Set founders in for-loop
SP <- SimParam

#Set Global Parameters
founderHap <- runMacs2(nInd=(nPa*4),nChr=10,Ne=Ne,segSites=1000)
SP <- SimParam$new(founderHap)
SP$restrSegSites(minQtlPerChr=1,minSnpPerChr=1,overlap=TRUE)
SP$addTraitA(nQtlPerChr=z1QTL,mean=z1mean,var=z1var)
SP$addTraitA(nQtlPerChr=z2QTL,mean=z2mean,var=z2var)
SP$setSexes("no")
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
SP$addSnpChip(nSnpPerChr=1000)
founders <- newPop(founderHap,simParam=SP)
founders <- setPheno(pop=founders,h2=h,reps=1)

## Set founding population breeding and selection parameters
parents <- nPa*2
crosses <- ceiling(factorial(parents)/(factorial(2)*factorial(parents-2)))
progeny <- nPr*2

## Set Breeding and Selection Parameters
nCycles <- 40 # loop length / number of cycles of selection
nParents <- ceiling(nPa) # number of parents to select and cross in first cycle
nCrosses <- ceiling(factorial(nParents)/(factorial(2)*factorial(nParents-2))) #nCr - total number of reciprocal crosses given nParents
nProgeny <- ceiling(nPr) # number of progeny taken from each parent

chosenParents <- selectInd(founders,nInd=parents,use="rand")
offspringPop <- randCross(pop=chosenParents,nCrosses=crosses,nProgeny=progeny)
offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)

simOutputWild <- list(founders,offspringPop)
simOutputOrphan <- list(founders,offspringPop)
simOutputLandrace <- list(founders,offspringPop)

rm(founders,founderHap,offspringPop)

for (pop in popType) { 
  
  
  for (cycle in 2:nCycles) { # Burn-in phase (creating populations)
    
    cat(paste0(" C",cycle))
    if (pop == "wild") chosenParents <- selectInd(simOutputWild[[cycle]],nInd=ceiling(length(simOutputWild[[cycle]]@id)/2),use="rand")
    if (pop == "wild") offspringPop <- randCross(pop=chosenParents,nCrosses=crosses,nProgeny=progeny)
    if (pop == "wild") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "wild") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan") allParents <- as.data.frame(simOutputOrphan[[cycle]]@pheno)
    if (pop == "orphan") allParents$id <- simOutputOrphan[[cycle]]@id
    if (pop == "orphan") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan") si <- (allParents[1]*0.5)+(allParents[2]*0.5)
    if (pop == "orphan") allParents$si <- si$Trait1
    if (pop == "orphan") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "orphan") chosenParents <- ordParents[1:(ceiling(length(simOutputOrphan[[cycle]]@id)*.75)),]
    if (pop == "orphan") chosenParents <- simOutputOrphan[[cycle]][chosenParents$id]
    if (pop == "orphan") offspringPop <- randCross(pop=chosenParents,nCrosses=crosses,nProgeny=progeny)
    if (pop == "orphan") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "orphan") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace") allParents <- as.data.frame(simOutputLandrace[[cycle]]@pheno)
    if (pop == "landrace") allParents$id <- simOutputLandrace[[cycle]]@id
    if (pop == "landrace") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace") si <- (allParents[1]*0.5)+(allParents[2]*0.5)
    if (pop == "landrace") allParents$si <- si$Trait1
    if (pop == "landrace") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "landrace") chosenParents <- ordParents[1:(ceiling(length(simOutputLandrace[[cycle]]@id)*.5)),]
    if (pop == "landrace") chosenParents <- simOutputLandrace[[cycle]][chosenParents$id]
    if (pop == "landrace") offspringPop <- randCross(pop=chosenParents,nCrosses=crosses,nProgeny=progeny)
    if (pop == "landrace") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "landrace") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
  }
  
  for (cycle in 41:42) {
    cat(paste0(" C",cycle))
    # PRS - Phenotypic Recurrent Selection (PRS following Germplasm Acquisition for 2 years)
    if (pop == "wild" & sel == "PRS") allParents <- as.data.frame(simOutputWild[[cycle]]@pheno)
    if (pop == "wild" & sel == "PRS") allParents$id <- simOutputWild[[cycle]]@id
    if (pop == "wild" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "wild" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "wild" & sel == "PRS") chosenParents <- simOutputWild[[cycle]][chosenParents$id]
    if (pop == "wild" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "wild" & sel == "PRS") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan" & sel == "PRS") allParents <- as.data.frame(simOutputOrphan[[cycle]]@pheno)
    if (pop == "orphan" & sel == "PRS") allParents$id <- simOutputOrphan[[cycle]]@id
    if (pop == "orphan" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "orphan" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "orphan" & sel == "PRS") chosenParents <- simOutputOrphan[[cycle]][chosenParents$id]
    if (pop == "orphan" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "orphan" & sel == "PRS") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace" & sel == "PRS") allParents <- as.data.frame(simOutputLandrace[[cycle]]@pheno)
    if (pop == "landrace" & sel == "PRS") allParents$id <- simOutputLandrace[[cycle]]@id
    if (pop == "landrace" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "landrace" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "landrace" & sel == "PRS") chosenParents <- simOutputLandrace[[cycle]][chosenParents$id]
    if (pop == "landrace" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "landrace" & sel == "PRS") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
    # GS - Genomic Selection (PRS following Germplasm Acquisition for 2 years)
    if (pop == "wild" & sel == "GS") allParents <- as.data.frame(simOutputWild[[cycle]]@pheno)
    if (pop == "wild" & sel == "GS") allParents$id <- simOutputWild[[cycle]]@id
    if (pop == "wild" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "GS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "wild" & sel == "GS") chosenParents <- ordParents[1:nParents,]
    if (pop == "wild" & sel == "GS") chosenParents <- simOutputWild[[cycle]][chosenParents$id]
    if (pop == "wild" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "GS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "wild" & sel == "GS") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan" & sel == "GS") allParents <- as.data.frame(simOutputOrphan[[cycle]]@pheno)
    if (pop == "orphan" & sel == "GS") allParents$id <- simOutputOrphan[[cycle]]@id
    if (pop == "orphan" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "GS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "orphan" & sel == "GS") chosenParents <- ordParents[1:nParents,]
    if (pop == "orphan" & sel == "GS") chosenParents <- simOutputOrphan[[cycle]][chosenParents$id]
    if (pop == "orphan" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "GS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "orphan" & sel == "GS") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace" & sel == "GS") allParents <- as.data.frame(simOutputLandrace[[cycle]]@pheno)
    if (pop == "landrace" & sel == "GS") allParents$id <- simOutputLandrace[[cycle]]@id
    if (pop == "landrace" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "GS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "landrace" & sel == "GS") chosenParents <- ordParents[1:nParents,]
    if (pop == "landrace" & sel == "GS") chosenParents <- simOutputLandrace[[cycle]][chosenParents$id]
    if (pop == "landrace" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "GS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "landrace" & sel == "GS") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
    # MxAv - Maximum Avoidance Selection (PRS following Germplasm Acquisition for 2 years)
    if (pop == "wild" & sel == "MxAv") allParents <- as.data.frame(simOutputWild[[cycle]]@pheno)
    if (pop == "wild" & sel == "MxAv") allParents$id <- simOutputWild[[cycle]]@id
    if (pop == "wild" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "MxAv") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "wild" & sel == "MxAv") chosenParents <- ordParents[1:nParents,]
    if (pop == "wild" & sel == "MxAv") chosenParents <- simOutputWild[[cycle]][chosenParents$id]
    if (pop == "wild" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "MxAv") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "wild" & sel == "MxAv") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan" & sel == "MxAv") allParents <- as.data.frame(simOutputOrphan[[cycle]]@pheno)
    if (pop == "orphan" & sel == "MxAv") allParents$id <- simOutputOrphan[[cycle]]@id
    if (pop == "orphan" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "MxAv") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "orphan" & sel == "MxAv") chosenParents <- ordParents[1:nParents,]
    if (pop == "orphan" & sel == "MxAv") chosenParents <- simOutputOrphan[[cycle]][chosenParents$id]
    if (pop == "orphan" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "MxAv") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "orphan" & sel == "MxAv") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace" & sel == "MxAv") allParents <- as.data.frame(simOutputLandrace[[cycle]]@pheno)
    if (pop == "landrace" & sel == "MxAv") allParents$id <- simOutputLandrace[[cycle]]@id
    if (pop == "landrace" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "MxAv") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "landrace" & sel == "MxAv") chosenParents <- ordParents[1:nParents,]
    if (pop == "landrace" & sel == "MxAv") chosenParents <- simOutputLandrace[[cycle]][chosenParents$id]
    if (pop == "landrace" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "MxAv") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "landrace" & sel == "MxAv") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
  }
  
  for (cycle in nCycles+3:nCycles) { # Selection phase (applying different selection techniques)
    cat(paste0(" C",cycle))
    # PRS - Phenotypic Recurrent Selection
    if (pop == "wild" & sel == "PRS") allParents <- as.data.frame(simOutputWild[[cycle]]@pheno)
    if (pop == "wild" & sel == "PRS") allParents$id <- simOutputWild[[cycle]]@id
    if (pop == "wild" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "wild" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "wild" & sel == "PRS") chosenParents <- simOutputWild[[cycle]][chosenParents$id]
    if (pop == "wild" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "wild" & sel == "PRS") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan" & sel == "PRS") allParents <- as.data.frame(simOutputOrphan[[cycle]]@pheno)
    if (pop == "orphan" & sel == "PRS") allParents$id <- simOutputOrphan[[cycle]]@id
    if (pop == "orphan" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "orphan" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "orphan" & sel == "PRS") chosenParents <- simOutputOrphan[[cycle]][chosenParents$id]
    if (pop == "orphan" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "orphan" & sel == "PRS") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace" & sel == "PRS") allParents <- as.data.frame(simOutputLandrace[[cycle]]@pheno)
    if (pop == "landrace" & sel == "PRS") allParents$id <- simOutputLandrace[[cycle]]@id
    if (pop == "landrace" & sel == "PRS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "PRS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "PRS") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "PRS") ordParents <- allParents[with(allParents,order(allParents$si,decreasing=TRUE)),]
    if (pop == "landrace" & sel == "PRS") chosenParents <- ordParents[1:nParents,]
    if (pop == "landrace" & sel == "PRS") chosenParents <- simOutputLandrace[[cycle]][chosenParents$id]
    if (pop == "landrace" & sel == "PRS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "PRS") offspringPop <- setPheno(pop=offspringPop,h2=h,reps=1)
    if (pop == "landrace" & sel == "PRS") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
    # GS - Genomic Selection (retraining every 2 years)
    nYearsToUse <- 2
    if (pop == "wild" & sel == "GS") trainingCyclesOutputWild <-purrr::reduce(simOutputWild %>% tail(nYearsToUse),`c`)
    if (pop == "wild" & sel == "GS") gen.mat <- pullSnpGeno(trainingCyclesOutputWild,simParam=SP)
    if (pop == "wild" & sel == "GS") A <- rrBLUP::A.mat(gen.mat)
    if (pop == "wild" & sel == "GS") allParents <- as.data.frame(trainingCyclesOutputWild@pheno)
    if (pop == "wild" & sel == "GS") allParents$id <- trainingCyclesOutputWild@id
    if (pop == "wild" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "GS") ans_gs <- mmer(si~1,random=~vs(id,Gu=A),rcov=~units,data=allParents)
    if (pop == "wild" & sel == "GS") gebv <- tibble(id=names(ans_gs$U$`u:id`$si),gebv=as.numeric(ans_gs$U$`u:id`$si))
    if (pop == "wild" & sel == "GS") chosenParents <- gebv %>% arrange(desc(gebv)) %>% slice(1:nParents)
    if (pop == "wild" & sel == "GS") chosenParents <- trainingCyclesOutputWild[chosenParents$id]
    if (pop == "wild" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "GS") simOutputWild[[length(simOutputWild)]]<-setPheno(pop = simOutputWild[[length(simOutputWild)]],h2=h,reps=1)
    if (pop == "wild" & sel == "GS") simOutputWild[[length(simOutputWild)+1]] <- offspringPop
    
    if (pop == "orphan" & sel == "GS") trainingCyclesOutputOrphan <-purrr::reduce(simOutputOrphan %>% tail(nYearsToUse),`c`)
    if (pop == "orphan" & sel == "GS") gen.mat <- pullSnpGeno(trainingCyclesOutputOrphan,simParam=SP)
    if (pop == "orphan" & sel == "GS") A <- rrBLUP::A.mat(gen.mat)
    if (pop == "orphan" & sel == "GS") allParents <- as.data.frame(trainingCyclesOutputOrphan@pheno)
    if (pop == "orphan" & sel == "GS") allParents$id <- trainingCyclesOutputOrphan@id
    if (pop == "orphan" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "GS") ans_gs <- mmer(si~1,random=~vs(id,Gu=A),rcov=~units,data=allParents)
    if (pop == "orphan" & sel == "GS") gebv <- tibble(id=names(ans_gs$U$`u:id`$si),gebv=as.numeric(ans_gs$U$`u:id`$si))
    if (pop == "orphan" & sel == "GS") chosenParents <- gebv %>% arrange(desc(gebv)) %>% slice(1:nParents)
    if (pop == "orphan" & sel == "GS") chosenParents <- trainingCyclesOutputOrphan[chosenParents$id]
    if (pop == "orphan" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "GS") simOutputOrphan[[length(simOutputOrphan)]]<-setPheno(pop = simOutputOrphan[[length(simOutputOrphan)]],h2=h,reps=1)
    if (pop == "orphan" & sel == "GS") simOutputOrphan[[length(simOutputOrphan)+1]] <- offspringPop
    
    if (pop == "landrace" & sel == "GS") trainingCyclesOutputLandrace <-purrr::reduce(simOutputLandrace %>% tail(nYearsToUse),`c`)
    if (pop == "landrace" & sel == "GS") gen.mat <- pullSnpGeno(trainingCyclesOutputLandrace)
    if (pop == "landrace" & sel == "GS") A <- rrBLUP::A.mat(gen.mat)
    if (pop == "landrace" & sel == "GS") allParents <- as.data.frame(trainingCyclesOutputLandrace@pheno)
    if (pop == "landrace" & sel == "GS") allParents$id <- trainingCyclesOutputLandrace@id
    if (pop == "landrace" & sel == "GS") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "GS") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "GS") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "GS") ans_gs <- mmer(si~1,random=~vsr(id,Gu=A),rcov=~units,data=allParents)
    if (pop == "landrace" & sel == "GS") gebv <- tibble(id=names(ans_gs$U$`u:id`$si),gebv=as.numeric(ans_gs$U$`u:id`$si))
    if (pop == "landrace" & sel == "GS") chosenParents <- gebv %>% arrange(desc(gebv)) %>% slice(1:nParents)
    if (pop == "landrace" & sel == "GS") chosenParents <- trainingCyclesOutputLandrace[chosenParents$id]
    if (pop == "landrace" & sel == "GS") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "GS") simOutputLandrace[[length(simOutputLandrace)]]<-setPheno(pop = simOutputLandrace[[length(simOutputLandrace)]],h2=h,reps=1)
    if (pop == "landrace" & sel == "GS") simOutputLandrace[[length(simOutputLandrace)+1]] <- offspringPop
    
    # MxAv - Maximum Avoidance Selection
    if (pop == "wild" & sel == "MxAv") trainingCyclesOutputWild <-purrr::reduce(simOutputWild %>% tail(nYearsToUse),`c`)
    if (pop == "wild" & sel == "MxAv") allParents <- as.data.frame(trainingCyclesOutputWild@pheno)
    if (pop == "wild" & sel == "MxAv") allParents$id <- trainingCyclesOutputWild@id
    if (pop == "wild" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "wild" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "wild" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "wild" & sel == "MxAv") id <- trainingCyclesOutputWild@id
    if (pop == "wild" & sel == "MxAv") dadid <- trainingCyclesOutputWild@father
    if (pop == "wild" & sel == "MxAv") momid <- trainingCyclesOutputWild@mother
    if (pop == "wild" & sel == "MxAv") pKin <- kinship2::kinship(id=id,dadid=dadid,momid=momid)*2
    if (pop == "wild" & sel == "MxAv") ans_MxAv <- mmer(si~1,random=~vs(id,Gu=pKin),rcov=~units,data=allParents)
    if (pop == "wild" & sel == "MxAv") pebv <- tibble(id=names(ans_MxAv$U$`u:id`$si),gebv=as.numeric(ans_MxAv$U$`u:id`$si))
    if (pop == "wild" & sel == "MxAv") chosenParents <- pebv %>% arrange(desc(gebv)) %>% slice(1:(nParents*4))
    if (pop == "wild" & sel == "MxAv") pKinuse <- pKin[rownames(pKin) %in% chosenParents$id,colnames(pKin) %in% chosenParents$id]
    if (pop == "wild" & sel == "MxAv") phen <- cbind(chosenParents$id,chosenParents$gebv)
    if (pop == "wild" & sel == "MxAv") colnames(phen) <- c("Indiv","BV")
    if (pop == "wild" & sel == "MxAv") phen<- as.data.frame(phen)
    if (pop == "wild" & sel == "MxAv") phen$Born <- c(cycle)
    if (pop == "wild" & sel == "MxAv") con  <- list()
    if (pop == "wild" & sel == "MxAv") cand <- candes(phen=phen, pKin=pKinuse, cont=NULL)
    if (pop == "wild" & sel == "MxAv") Offspring <- opticont("min.pKin", cand, con,solver="alabama")
    if (pop == "wild" & sel == "MxAv") out <- Offspring$parent$Indiv
    if (pop == "wild" & sel == "MxAv") MxAvSel <- out[1:nParents]
    if (pop == "wild" & sel == "MxAv") chosenParents <- trainingCyclesOutputWild[MxAvSel]
    if (pop == "wild" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "wild" & sel == "MxAv") simOutputWild[[length(simOutputWild)]]<-setPheno(pop = simOutputWild[[length(simOutputWild)]],h2=h,reps=1)
    if (pop == "wild" & sel == "MxAv") simOutputWild[[length(simOutputWild)+1]]<-offspringPop
    
    if (pop == "orphan" & sel == "MxAv") trainingCyclesOutputOrphan <-purrr::reduce(simOutputOrphan %>% tail(nYearsToUse),`c`)
    if (pop == "orphan" & sel == "MxAv") allParents <- as.data.frame(trainingCyclesOutputOrphan@pheno)
    if (pop == "orphan" & sel == "MxAv") allParents$id <- trainingCyclesOutputOrphan@id
    if (pop == "orphan" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "orphan" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "orphan" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "orphan" & sel == "MxAv") id <- trainingCyclesOutputOrphan@id
    if (pop == "orphan" & sel == "MxAv") dadid <- trainingCyclesOutputOrphan@father
    if (pop == "orphan" & sel == "MxAv") momid <- trainingCyclesOutputOrphan@mother
    if (pop == "orphan" & sel == "MxAv") pKin <- kinship2::kinship(id=id,dadid=dadid,momid=momid)*2
    if (pop == "orphan" & sel == "MxAv") ans_MxAv <- mmer(si~1,random=~vs(id,Gu=pKin),rcov=~units,data=allParents)
    if (pop == "orphan" & sel == "MxAv") pebv <- tibble(id=names(ans_MxAv$U$`u:id`$si),gebv=as.numeric(ans_MxAv$U$`u:id`$si))
    if (pop == "orphan" & sel == "MxAv") chosenParents <- pebv %>% arrange(desc(gebv)) %>% slice(1:(nParents*4))
    if (pop == "orphan" & sel == "MxAv") pKinuse <- pKin[rownames(pKin) %in% chosenParents$id,colnames(pKin) %in% chosenParents$id]
    if (pop == "orphan" & sel == "MxAv") phen <- cbind(chosenParents$id,chosenParents$gebv)
    if (pop == "orphan" & sel == "MxAv") colnames(phen) <- c("Indiv","BV")
    if (pop == "orphan" & sel == "MxAv") phen<- as.data.frame(phen)
    if (pop == "orphan" & sel == "MxAv") phen$Born <- c(cycle)
    if (pop == "orphan" & sel == "MxAv") con  <- list()
    if (pop == "orphan" & sel == "MxAv") cand <- candes(phen=phen, pKin=pKinuse, cont=NULL)
    if (pop == "orphan" & sel == "MxAv") Offspring <- opticont("min.pKin", cand, con,solver="alabama")
    if (pop == "orphan" & sel == "MxAv") out <- Offspring$parent$Indiv
    if (pop == "orphan" & sel == "MxAv") MxAvSel <- out[1:nParents]
    if (pop == "orphan" & sel == "MxAv") chosenParents <- trainingCyclesOutputOrphan[MxAvSel]
    if (pop == "orphan" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "orphan" & sel == "MxAv") simOutputOrphan[[length(simOutputOrphan)]]<-setPheno(pop = simOutputOrphan[[length(simOutputOrphan)]],h2=h,reps=1)
    if (pop == "orphan" & sel == "MxAv") simOutputOrphan[[length(simOutputOrphan)+1]]<-offspringPop
    
    if (pop == "landrace" & sel == "MxAv") trainingCyclesOutputLandrace <-purrr::reduce(simOutputLandrace %>% tail(nYearsToUse),`c`)
    if (pop == "landrace" & sel == "MxAv") allParents <- as.data.frame(trainingCyclesOutputLandrace@pheno)
    if (pop == "landrace" & sel == "MxAv") allParents$id <- trainingCyclesOutputLandrace@id
    if (pop == "landrace" & sel == "MxAv") allParents[,1:2] <- scale(allParents[,1:2])
    if (pop == "landrace" & sel == "MxAv") si <- (allParents[1]*siC)+(allParents[2]*(1-siC))
    if (pop == "landrace" & sel == "MxAv") allParents$si <- si$Trait1
    if (pop == "landrace" & sel == "MxAv") id <- trainingCyclesOutputLandrace@id
    if (pop == "landrace" & sel == "MxAv") dadid <- trainingCyclesOutputLandrace@father
    if (pop == "landrace" & sel == "MxAv") momid <- trainingCyclesOutputLandrace@mother
    if (pop == "landrace" & sel == "MxAv") pKin <- kinship2::kinship(id=id,dadid=dadid,momid=momid)*2
    if (pop == "landrace" & sel == "MxAv") ans_MxAv <- mmer(si~1,random=~vs(id,Gu=pKin),rcov=~units,data=allParents)
    if (pop == "landrace" & sel == "MxAv") pebv <- tibble(id=names(ans_MxAv$U$`u:id`$si),gebv=as.numeric(ans_MxAv$U$`u:id`$si))
    if (pop == "landrace" & sel == "MxAv") chosenParents <- pebv %>% arrange(desc(gebv)) %>% slice(1:(nParents*4))
    if (pop == "landrace" & sel == "MxAv") pKinuse <- pKin[rownames(pKin) %in% chosenParents$id,colnames(pKin) %in% chosenParents$id]
    if (pop == "landrace" & sel == "MxAv") phen <- cbind(chosenParents$id,chosenParents$gebv)
    if (pop == "landrace" & sel == "MxAv") colnames(phen) <- c("Indiv","BV")
    if (pop == "landrace" & sel == "MxAv") phen<- as.data.frame(phen)
    if (pop == "landrace" & sel == "MxAv") phen$Born <- c(cycle)
    if (pop == "landrace" & sel == "MxAv") con  <- list()
    if (pop == "landrace" & sel == "MxAv") cand <- candes(phen=phen, pKin=pKinuse, cont=NULL)
    if (pop == "landrace" & sel == "MxAv") Offspring <- opticont("min.pKin", cand, con,solver="alabama")
    if (pop == "landrace" & sel == "MxAv") out <- Offspring$parent$Indiv
    if (pop == "landrace" & sel == "MxAv") MxAvSel <- out[1:nParents]
    if (pop == "landrace" & sel == "MxAv") chosenParents <- trainingCyclesOutputLandrace[MxAvSel]
    if (pop == "landrace" & sel == "MxAv") offspringPop <- randCross(pop=chosenParents,nCrosses=nCrosses,nProgeny=nProgeny)
    if (pop == "landrace" & sel == "MxAv") simOutputLandrace[[length(simOutputLandrace)]]<-setPheno(pop = simOutputLandrace[[length(simOutputLandrace)]],h2=h,reps=1)
    if (pop == "landrace" & sel == "MxAv") simOutputLandrace[[length(simOutputLandrace)+1]]<-offspringPop
    
  }
  if (pop == "wild") tidySimOutputWild<-tibble(Cycle=0:(length(simOutputWild)-1)) %>% 
      dplyr::mutate(meanPz1=map_dbl(simOutputWild,~mean(.@pheno[,1])),
                    varPz1=map_dbl(simOutputWild,~var(.@pheno[,1])),
                    meanPz2=map_dbl(simOutputWild,~mean(.@pheno[,2])),
                    varPz2=map_dbl(simOutputWild,~var(.@pheno[,2])),
                    Nprogeny=nProgeny,
                    Nparents=nParents,
                    NeStart=Ne)
  if (pop == "wild") tidySimOutputWild$si<-rep(paste(siC,collapse=''),nrow(tidySimOutputWild))
  if (pop == "wild") tidySimOutputWild$h2<-rep(paste(h,collapse=''),nrow(tidySimOutputWild))
  if (pop == "wild") tidySimOutputWild$sel<-rep(paste(sel,collapse=''),nrow(tidySimOutputWild))
  if (pop == "wild") tidySimOutputWild$pop<-rep(paste(pop,collapse=''),nrow(tidySimOutputWild))
  
  if (pop == "orphan") tidySimOutputOrphan<-tibble(Cycle=0:(length(simOutputOrphan)-1)) %>% 
      dplyr::mutate(meanPz1=map_dbl(simOutputOrphan,~mean(.@pheno[,1])),
                    varPz1=map_dbl(simOutputOrphan,~var(.@pheno[,1])),
                    meanPz2=map_dbl(simOutputOrphan,~mean(.@pheno[,2])),
                    varPz2=map_dbl(simOutputOrphan,~var(.@pheno[,2])),
                    Nprogeny=nProgeny,
                    Nparents=nParents,
                    NeStart=Ne)
  if (pop == "orphan") tidySimOutputOrphan$si<-rep(paste(siC,collapse=''),nrow(tidySimOutputOrphan))
  if (pop == "orphan") tidySimOutputOrphan$h2<-rep(paste(h,collapse=''),nrow(tidySimOutputOrphan))
  if (pop == "orphan") tidySimOutputOrphan$sel<-rep(paste(sel,collapse=''),nrow(tidySimOutputOrphan))
  if (pop == "orphan") tidySimOutputOrphan$pop<-rep(paste(pop,collapse=''),nrow(tidySimOutputOrphan))
  
  if (pop == "landrace") tidySimOutputLandrace<-tibble(Cycle=0:(length(simOutputLandrace)-1)) %>% 
      dplyr::mutate(meanPz1=map_dbl(simOutputLandrace,~mean(.@pheno[,1])),
                    varPz1=map_dbl(simOutputLandrace,~var(.@pheno[,1])),
                    meanPz2=map_dbl(simOutputLandrace,~mean(.@pheno[,2])),
                    varPz2=map_dbl(simOutputLandrace,~var(.@pheno[,2])),
                    Nprogeny=nProgeny,
                    Nparents=nParents,
                    NeStart=Ne)
  if (pop == "landrace") tidySimOutputLandrace$si<-rep(paste(siC,collapse=''),nrow(tidySimOutputLandrace))
  if (pop == "landrace") tidySimOutputLandrace$h2<-rep(paste(h,collapse=''),nrow(tidySimOutputLandrace))
  if (pop == "landrace") tidySimOutputLandrace$sel<-rep(paste(sel,collapse=''),nrow(tidySimOutputLandrace))
  if (pop == "landrace") tidySimOutputLandrace$pop<-rep(paste(pop,collapse=''),nrow(tidySimOutputLandrace))
}


# Combine Output lists to dataframes
all <- data.frame(rbind(tidySimOutputWild,tidySimOutputOrphan,tidySimOutputLandrace)) # Bind all population replications to 1 dataframe
all$h2 <- as.factor(all$h2)
all <- all %>% mutate(h2=recode(h2,"0.30.3"="hLL","0.30.7"="hLH","0.50.5"="hMM","0.70.3"="hHL")) # rename heritability used

# Get the SLURM_ARRAY_TASK_ID from the environment
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Construct the filename based on the task ID
output_file <- paste0("output_", task_id, ".csv")

# Save the data to a CSV file with the task ID in the filename
write.csv(all, file = output_file)



