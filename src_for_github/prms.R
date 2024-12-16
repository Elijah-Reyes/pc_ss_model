
make.prms <- function(Nm.init = 100,
                      Nf.init = 100,
                      male.care.inits = c(9, 0),
                      female.care.inits = c(0, 9),
                      male.sex.inits =  c(1,7),
                      female.sex.inits =  c(7,1),
                      L = 10,
                      num.traits = 5,
                      mu = 1e-4,
                      birth.ratio = 0.5,
                      lek = 1,
                      k = c(200, 200),
                      r = 10,
                      c = 0.75,
                      a = 0.05,
                      b = 1.5,
                      rho.care = 1.2,
                      tau = 0.15,
                      brood.size = 1,
                      care.saturate = 0.6,
                      rest.days = 0, 
                      mature.time = 10,
                      n.patches = 1,
                      m.disp = 0,
                      f.disp = 0,
                      mate.system = 'CTM',
                      mort.state=FALSE,
                      ss.off=FALSE,
                      store.fit = FALSE,
                      store.growth=FALSE,
                      store.hist=FALSE,
                      inits.fixed=TRUE,
                      fitness.ngens = 1000,
                      num.gens = 10000,
                      run.watch = FALSE,
                      print.every=100,
                      save.every=100,
                      run.index=NULL
                      ){
  
    inputs <- as.list(environment())
    inputs <- inputs[order(names(inputs))] ## order for consistency

    return(inputs)
}

## ~~~~~~~~~~~~~~~~~~~~
## ~~ Parameters ~~~~~~
## ~~~~~~~~~~~~~~~~~~~~
##
## Here is a summary of all present paramters and what they do:
##
## ~~ Initial pop ~~
##
## Nm.init:          initial size of the male population
## Nf.init:          initial size of the female population
## male.care.init:   initial proportion of care loci that are a '1' (male)
## female.care.init: initial proportion of care loci that are a '1' (female) 
## male.sex.init:    initial proportion of sexual loci that are a '1' (male)
## female.sex.init:  initial proportion of sexual loci that are a '1' (female)
## L:                number of loci for each trait
##
##
## ~~ Function params ~~
##
## mu:            mutation rate
## birth.ratio:   sex ratio of juveniles who survived. Actually MSR but
##                  juvenile survival is not sex biased so it is the same 
## lek:           number of males a female samples if females are choosy 
## k:             carrying capacity of the patch
## r:             growth rate
## c:             offspring survival coefficient
## brood.size:    mean number of offspring produced from each mating
##                  if pop is not density dependent 
## a:             number of timesteps to reenter mating pool outside
##                  of parental care time.
## mature.time:   number of timesteps for offspring to mature 
## mate.symmetry: male and female traits represent the same thing
##
##
## ~~ Running params ~~
##
## num.gens:  the number of generations to run the model
## run.watch: do you want to print and watch the population stats and
##            plots as it runs. Turn off if running parallel. 


