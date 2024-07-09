
## ~~~~~~~~~~~~~~~~~~~~~
## ~~ Source file ~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~
##
## Sources the other files that will be used in the model run. These
## are all found in the 'src' folder. These each do the following:
## prms:      Parameters function
## pop_fxns:  Functions for thing the lifecylce uses. For example, it
##              contains the reprodcution and death functions.
## lifecycle: Calls the functions for the population. Sets the order
##              of events that individuals experience.
## meta_fxns: A function that makes all of the intial objects and then
##              loops through generations. 

source('src/prms.R')
source('src/pop_fxns.R')
source('src/lifecycle.R')
source('src/meta_fxns.R')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Population function ~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## Function that makes the population. Resulting populaption is a list
## with 'males', 'females', and then an 'id' tracker that is used to
## keep track of individuals.  

make.pop <- function(prms) {
 
    males <- array(0, dim = c(prms$Nm.init, (prms$L * prms$num.traits) + 4), 
                   dimnames = list(1:prms$Nm.init, 
                                   c(sprintf('locus.%s', 1:(prms$L * prms$num.traits)), 
                                   'sex', 'ID', 'state', 'patch')))

    females <- array(0, dim = c(prms$Nf.init, (prms$L * prms$num.traits) + 4), 
                     dimnames = list(1:prms$Nf.init, 
                                     c(sprintf('locus.%s', 1:(prms$L * prms$num.traits)), 
                                     'sex', 'ID', 'state', 'patch')))


    fill.pop <- function(sex, traits, N.init, sex.num, prms=prms){
        for(ll in 1:prms$num.traits){
            sex[,(prms$L * (ll-1) + 1):(prms$L * ll)] <- rbinom(prms$L *
                                                                  N.init, 1,
                                                                  prob=traits[1,ll]) 
        }    
  
        sex[,'state'] <- 1
        sex[,'sex'] <- sex.num
        sex[,'patch'] <- 1
        sex[,'ID'] <- 1:N.init

        if(prms$n.patches==2){
            sex.2 <- sex
            sex.2[,'patch'] <- 2
            sex.2[,'ID'] <- (N.init + 1):(N.init*2)

            for(ll in 1:prms$num.traits){
                sex.2[,(prms$L * (ll-1) + 1):(prms$L * ll)] <- rbinom(prms$L *
                                                                        N.init, 1,
                                                                        prob=traits[2,ll]) 
            }

            sex <- rbind(sex, sex.2)
        }

        pheno.sexm <- rowSums(sex[,1:prms$L]) 
        pheno.pcm  <- rowSums(sex[,(prms$L + 1):(prms$L * 2)])
        pheno.sexf <- rowSums(sex[,(prms$L * 2 + 1):(prms$L * 3)]) 
        pheno.pcf  <- rowSums(sex[,(prms$L * 3 + 1):(prms$L * 4)])

        sex <- cbind(sex,
                     pheno.sexm, pheno.pcm,
                     pheno.sexf, pheno.pcf)

        return(sex)
    }

    
    traits <- rbind(c(prms$male.sex.inits[1], prms$male.care.inits[1], prms$female.sex.inits[1], prms$female.care.inits[1], 0), 
                    c(prms$male.sex.inits[2], prms$male.care.inits[2], prms$female.sex.inits[2], prms$female.care.inits[2], 1))


    males <- fill.pop(sex=males, traits=traits, N.init=prms$Nm.init, sex.num=1, prms=prms)
    females <- fill.pop(sex=females, traits=traits, N.init=prms$Nf.init, sex.num=0, prms=prms)
  
    pop <- list(females = females,
                males = males,
                id.counter = c(male.id   = nrow(males),
                               female.id = nrow(females),
                               juvenile.id = 0))
    return(pop)
}

## Make an objects for storing the immature juveniles
make.juv <- function(){
    juveniles <- array(NA, dim=c(0,8))
    dimnames(juveniles)[[2]] <- c('mom.id', 'dad.id',
                                  'mom.ind','dad.ind',
                                  'mom.care','dad.care',
                                  'mature', 'patch')
    return(juveniles)
}

## make an object to store the adults who have died and have living immature offspring
make.coffin <- function(pp=pop$males){
    dead.f <- array(NA, dim=c(0,ncol(pp))) 
    dimnames(dead.f)[[2]] <- dimnames(pp)[[2]]
    dead.m <- array(NA, dim=c(0,ncol(pp)))
    dimnames(dead.m)[[2]] <- dimnames(pp)[[2]]

    coffin <- list(dead.f = dead.f, dead.m = dead.m)
    return(coffin)
}

## make an empty object to store fitness measures 
make.fitbit <- function(pop=pop, prms=prms){
    fit.m <- array(NA, dim=c(0,5))
    dimnames(fit.m) <- list(NULL, c('sex', 'ID', 'patch', 'parCare',
                             'sexSel'))
    fit.f <- fit.m

    fit.stats <- array(NA, dim=c(2,2,2,3), 
                    dimnames=list(sex=c('male', 'female'),
                                patch=c(1,2),
                                patch.natal=c(1,2),
                                fitness=c('matings', 'off.survival','total.offspring')))                               

    fitbit <- list(fit.m = fit.m, fit.f = fit.f, fit.stats=fit.stats)

    fitbit$fit.m <- rbind(fitbit$fit.m,
                        pop$males[,c('sex', 'ID',
                                     'patch',
                                     'pheno.pcm',
                                     'pheno.sexm')])
    fitbit$fit.m <- cbind(fitbit$fit.m, matings=0, off.surv=0,
                        patch.natal=0, total.off=0, lifespan=0) 

    fitbit$fit.f <- rbind(fitbit$fit.f,
                        pop$females[,c('sex', 'ID',
                                       'patch',
                                       'pheno.pcf',
                                       'pheno.sexf')])
    fitbit$fit.f <- cbind(fitbit$fit.f, matings=0, off.surv=0,
                        patch.natal=0, total.off=0, lifespan=0)
    
    return(fitbit)
}