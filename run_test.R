## Set working directory here:
setwd()

## Resets environment
rm(list=ls())

## Sources model
source('src/init.R')

## Runs model; see prms.R for default parameters; 
## any parameters in the prms.R file can be updated here as an arguement to run.model()    

out <- run.model(male.care.inits = c(0.9,0.05),   
                 female.care.inits = c(0.05,0.9),   
                 male.sex.inits = c(0.05,0.9),    
                 female.sex.inits = c(0.9,0.05),
                 n.patches = 2,  
                 mate.system='limited.sex', ##Choose 1: 'limited.sex' 'limited.choosy', choosy.sex
                 num.gens=1e4,
                 store.fit=TRUE,
                 store.growth=TRUE,
                 fitness.ngens=1000,
                 ss.off=F,
                 run.watch=TRUE,
                 print.every=1) 
