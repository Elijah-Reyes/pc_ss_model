
lifecycle <- function(pop, juveniles, coffin, prms, fitbit, pop.growth, gen){

## figures out who survives 
 if(prms$mort.state==TRUE){
    f.surv <- mortality.state(pp=pop$females,
                              prms=prms)
    m.surv <- mortality.state(pp=pop$males,
                              prms=prms)   
 }else{
    f.surv <- mortality(pp=pop$females,
                        prms=prms)
    m.surv <- mortality(pp=pop$males,
                        prms=prms)
 }
                 
## storing adults that have died that have offspring maturing
    coffin$dead.f <- update.coffin(surv=f.surv, foc.pop=pop$females, foc.dead=coffin$dead.f, 
                                 juvs=juveniles, foc.pop.id='mom.id')
    coffin$dead.m <- update.coffin(surv=m.surv, foc.pop=pop$males, foc.dead=coffin$dead.m, 
                                   juvs=juveniles, foc.pop.id='dad.id')
    pop$females <- pop$females[f.surv == 1,,drop=FALSE]
    pop$males <- pop$males[m.surv == 1,,drop=FALSE]
   
    if(prms$store.growth==TRUE){
        pop.growth['adult.deaths',gen] <- sum(f.surv == 0, m.surv == 0)
    }

## figures out who mates and produces their offspring. Goes through all patches.
    parents <- vector()

    for(pp in 1:prms$n.patches){

        foc.mm <- pop$males[,'patch']==pp
        foc.ff <- pop$females[,'patch']==pp
    
        foc.par <- mating(mm=pop$males,
                          ff=pop$females,
                          prms,
                          patch=pp)

        pop$females[foc.par[,'mom.ind'],'state'] <- 2
        pop$males[foc.par[,'dad.ind'],'state'] <- 2

        foc.par <- cbind(foc.par,
                    mom.care = pop$females[foc.par[,'mom.ind'],'pheno.pcf'],
                    dad.care =
                        pop$males[foc.par[,'dad.ind'],'pheno.pcm'],
                    mature = rep(prms$mature.time, nrow(foc.par)),
                    patch=rep(pp, nrow(foc.par)))

        if(nrow(foc.par) > 0){
        broods <- rpois(nrow(foc.par),
                        lambda = max(fecundity(N = (sum(foc.mm) +
                                                    sum(foc.ff)),
                                            prms=prms,
                                            patch=pp), 0)) 

        if(prms$store.fit==TRUE & gen < (prms$fitness.ngens+1)){
            fitbit <- update.fitbit(fitbit=fitbit, parents=foc.par, broods=broods)
        } 

        foc.par <- foc.par[rep(1:nrow(foc.par), broods),,drop=FALSE]
        }
        
        parents <- rbind(parents, foc.par)
    }
                                                                  
## stiches juveniles and new offspring parent data together for passing across timesteps
    juveniles <- rbind(juveniles, parents)
    if(prms$store.growth==TRUE){
        pop.growth['births',gen] <- nrow(parents)
    }
    off.surv <- off.survival(juveniles, prms)
    survived <- rbinom(length(off.surv), 1, prob=off.surv)
    juveniles <- juveniles[survived > 0,,drop=FALSE]
    if(prms$store.growth==TRUE){
        pop.growth['offspring.deaths',gen] <- length(off.surv) - sum(survived)
    }

## update parents without offspring to move to mating pool  
    pop$females[(pop$females[,'ID'] %in% juveniles[,'mom.id']) == FALSE,
                'state'] <- 1
    pop$males[(pop$males[,'ID'] %in% juveniles[,'dad.id']) == FALSE,
                'state'] <- 1

    pop$females[pop$females[,'ID'] %in%
                juveniles[juveniles[,'mom.care']==(prms$rest.days * -1), 'mom.id'], 'state'] <- 1
    pop$males[pop$males[,'ID'] %in%
                juveniles[juveniles[,'dad.care']==(prms$rest.days * -1), 'dad.id'], 'state'] <- 1

## A counter for parental care and maturation
   juveniles[,c('mom.care', 'dad.care', 'mature')] <-
    juveniles[,c('mom.care', 'dad.care', 'mature')] - 1

## Matured juveniles are moved to the mating pool
  matured.ind <- which(juveniles[,'mature'] == 0)
  if(prms$store.growth==TRUE){
    pop.growth['offspring.matured',gen] <- length(matured.ind)
  }
  if(length(matured.ind) > 0){
    matured <- juveniles[matured.ind,,drop=FALSE]
    updates <- update.adults(pop=pop, prms=prms,
                  coffin=coffin, matured=matured, gen)
    
    pop <- updates$pop

    if(prms$store.fit == TRUE  & gen < (prms$fitness.ngens+1)){
        par.f <- table(matured[,'mom.id'])
        par.m <- table(matured[,'dad.id'])
        f.inds <- as.numeric(names(par.f))
        m.inds <- as.numeric(names(par.m))
        idm.inds <- match(m.inds, fitbit$fit.m[,'ID'])
        idf.inds <- match(f.inds, fitbit$fit.f[,'ID'])
      
        m.inds <- !is.na(idm.inds)
        f.inds <- !is.na(idf.inds)

        fitbit$fit.f[idf.inds[f.inds],'off.surv'] <- fitbit$fit.f[idf.inds[f.inds],'off.surv'] + par.f[f.inds]
        fitbit$fit.m[idm.inds[m.inds],'off.surv'] <- fitbit$fit.m[idm.inds[m.inds],'off.surv'] + par.m[m.inds]
        
        fitbit$fit.m <- rbind(fitbit$fit.m, updates$new.males)
        fitbit$fit.f <- rbind(fitbit$fit.f, updates$new.females)

        fitbit <- update.lifespan(pop, fitbit)
    }
    
    juveniles <- juveniles[-matured.ind,,drop=FALSE]

  }
  
  return(list(pop=pop, juveniles=juveniles, coffin=coffin, fitbit=fitbit, pop.growth=pop.growth))
}






