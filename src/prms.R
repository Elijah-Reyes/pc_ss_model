
make.prms <- function(Nm.init = 100, ## Initial male pop size                  
                      Nf.init = 100, ## Initial female pop size
                      male.care.inits = c(0.9, 0.05),    ## initial male parental care in patch 1 and 2
                      female.care.inits = c(0.05, 0.9),  ## initial female parental care in patch 1 and 2
                      male.sex.inits =  c(0.1,0.7),      ## initial male sexual trait in patch 1 and 2
                      female.sex.inits =  c(0.7,0.1),    ## initial female sexual trait in patch 1 and 2
                      L = 10,             ## number of loci for each trait
                      num.traits = 5,     ## How many traits each individual carries
                      mu = 1e-4,          ## mutation rate
                      birth.ratio = 0.5,  ## birth sex ratio
                      k = c(200, 200),    ## Carrying capacity in patch 1 and 2
                      r = 10,             ## baseline number of offspring per mating
                      c = 0.75,           ## baseline offspring survival
                      a = 0.05,           ## baseline adult mortality
                      b = 1.5,            ## shape coefficent of adult mortality function
                      rho.care = 1.2,     ## additive/subadditive effect of two parent care (<2 = subadditive)
                      tau = 0.15,         ## benefit of uniparental care 
                      rest.days = 0,      ## lag time between returning to mating pool after care
                      mature.time = 10,   ## how many days for offspring to mature
                      n.patches = 1,      ## 1 = single patch model; 2 = two patch model (with dispersal)
                      m.disp = 0.005,     ## male dispersal probability
                      f.disp = 0.005,     ## female dispersal probability
                      mate.system = 'limited.sex',  ## mating system: CTM = 'limited.sex'; MPM = 'limited.choosy'; PDM = 'choosy.sex'
                      mort.state=FALSE,   ## state-specific mortality off (FALSE) or on (TRUE)
                      ss.off=FALSE,       ## Is the sexual trait neutral? (FALSE means the sexual trait increases mating)
                      store.fit = FALSE,  ## Store fitness measures in output?
                      store.growth=TRUE,  ## Store population growth meansures in output?
                      fitness.ngens = 1000,  ## How many generations to store of the fitness measures
                      num.gens = 10000,      ## How many generations to run the model
                      run.watch = FALSE,     ## Do you want to view plots of the output as the model runs?
                      print.every=100,       ## How often to print the genertion number and summary stats
                      run.index=NULL         ## Stores whatever index I'm using to run replicates in the parameter list to find in the output for verifying the runs later
                      ){
  
    inputs <- as.list(environment())
    inputs <- inputs[order(names(inputs))] ## order for consistency

    return(inputs)
}



