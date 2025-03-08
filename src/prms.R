
make.prms <- function(Nm.init = 100,                    ## number of initial males
                      Nf.init = 100,                    ## number of initial females
                      male.care.inits = c(9, 0),        ## initial male parental care in (patch1, patch2)
                      female.care.inits = c(0, 9),      ## initial female parental care in (patch1, patch2)
                      male.sex.inits =  c(1,7),         ## initial male sexual trait in (patch1, patch2)
                      female.sex.inits =  c(7,1),       ## initial female sexual trait in (patch1, patch2)
                      L = 10,                           ## loci per trait
                      num.traits = 5,                   ## number of traits to track
                      mu = 1e-4,                        ## mutation rate
                      birth.ratio = 0.5,                ## birth sex ratio                     
                      k = c(200, 200),                  ## carrying capcity of (patch1, patch2)
                      r = 10,                           ## growth rate
                      c = 0.75,                         ## baseline offspring survival
                      a = 0.05,                         ## scales mortality function
                      b = 1.5,                          ## exponent of mortality function
                      rho.care = 1.2,                   ## additive effect of both parents providing care
                      tau = 0.15,                       ## contribution to survival of a parent providing care
                      brood.size = 1,                   ## brood size
                      rest.days = 0,                    ## timesteps to reenter mating pool after providing care
                      mature.time = 10,                 ## number of timesteps for offspring to mature
                      n.patches = 1,                    ## how many patches (1 or 2)
                      m.disp = 0,                       ## male dispersal rate
                      f.disp = 0,                       ## female dispersal rate
                      mate.system = 'CTM',              ## mating system: noSS.CTM, noSS.MPM, 'CTM', 'MPM' (noSS means no sexual selection)
                      ss.off=FALSE,                     ## is sexual selection off?
                      store.fit = FALSE,                ## store the fitness metrics in output?
                      store.growth=FALSE,               ## store growth metrics in output?
                      store.hist=FALSE,                 ## store data for histgrams in output?
                      inits.fixed=TRUE,                 ## are the inits fixed in each population or is there variance?
                      fitness.ngens = 1000,             ## how many timesteps to store fitness metrics?
                      num.gens = 10000,                 ## number of timesteps to run model
                      run.watch = FALSE,                ## do you want to watch the model run with plots? (only for single runs)
                      print.every=100,                  ## how often to see messages about run progress in sinlge runs
                      save.every=100,                   ## how often to save population data
                      run.index=NULL                    ## place to store information about all combinations of parameters being run across cores (Note: I run multiple model by indexing an array containing all paramter combinations to run.)
                      ){
  
    inputs <- as.list(environment())
    inputs <- inputs[order(names(inputs))] ## order for consistency

    return(inputs)
}