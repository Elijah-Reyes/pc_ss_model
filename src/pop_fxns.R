
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Offspring survival ~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to transform parental care into offspring survival. Works
## by taking how many parents are caring in the current time step. We
## than assign an amount of survival based upon whether there is 0, 1,
## or 2 parents caring. 'c' acts to give a baseline amount of
## survival. We can think of it as telling how much survival depends
## on parental care vs the baseline. 
## To get FJ results, this needs diminishing returns such that 2
## parents providing care contributes less than that of 1 + 1 (e.g., 1
## parent provides 0.5 total while 2 parents provide 0.75 total).  

off.survival <- function(juveniles, prms){

  n.care <- ((juveniles[,'mom.care'] > 0) +
            (juveniles[,'dad.care'] > 0)) + 1
  
  rho <- c(0, 1, prms$rho.care)
  
  surv <- prms$c + prms$tau * rho[n.care]
  surv <- pmin(surv, 1)
  return(surv)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Fecundity ~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function for density dependent fecundity. Limits the population
## size with 'k' carrying capacity with growth rate 'r'. Density is
## based only on total adult population 

fecundity <- function(N, prms, patch){
  R <- prms$r * (1 - N/prms$k[patch])
  return(R)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Mortality ~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to figure out which adults die. Following FJ, adults in
## any timestep can die at a probability depending on their sexual
## trait. 'pp' is the pop being used. Since males and females are in
## different objects, this is sex specific. 

mortality <- function(pp, prms){

  if(nrow(pp) < 1){return(integer(0))}
  
  sex <- pp[1, 'sex', drop=FALSE]
  
  prob.arr <- array(NA, dim=c(1,nrow(pp)))
  
  if(sex == 0){
    sex.phenos <- c('pheno.sexf','pheno.pcf')
  }else{
    sex.phenos <- c('pheno.sexm','pheno.pcm')
  }

  baseline.trait.mort <- prms$a * (1 + (pp[,sex.phenos[1]]/prms$L)^(prms$b)) 
  surv.prob <- 1 - baseline.trait.mort 
   
  surv.ind <- rbinom(nrow(pp), 1, prob=surv.prob)
  
  return(surv.ind)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Mating ~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Figures out which males and females are mating in the mating
## pool. See comments throughout for greater detail. 

alpha <- function(sex.trait, s){
    prob.f.mates <- (sex.trait*s)/10
    pmin(1, prob.f.mates)
}

mating <- function(mm, ff, prms, patch){

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##     ~~ Limiting sex ~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Figure out which sex is limiting. 

  m.ind <- which(mm[,'state'] == 1 & mm[,'patch'] == patch)
  f.ind <- which(ff[,'state'] == 1 & ff[,'patch'] == patch)
  if(length(m.ind)>1){m.ind <- sample(m.ind)}  
  if(length(f.ind)>1){f.ind <- sample(f.ind)}

  names(f.ind) <- names(m.ind) <- NULL
  
  if(length(m.ind)<1 | length(f.ind)<1){
    tmp <- array(NA,dim=c(0,4))
    dimnames(tmp)[[2]] <- c('mom.id', 'dad.id', 'mom.ind', 'dad.ind')
    return(tmp)
  }
  male.biased <- length(f.ind) < length(m.ind)

    if(prms$mate.system == 'noSS.CTM'){
        n.matings <- pmin(length(f.ind), length(m.ind))
        mated <- vector(length = n.matings)
        if(length(f.ind)==1 & length(m.ind)==1){
            moms <- f.ind
            dads <- m.ind
        }else{
            if(male.biased){
                avail.mates <- m.ind
                mated[] <- sample(avail.mates, n.matings,
                                    replace=FALSE)
                moms <- f.ind
                dads <- mated
            }else{
                avail.mates <- f.ind
                mated[] <- sample(avail.mates, n.matings,
                                    replace=FALSE) 
                moms <- mated
                dads <- m.ind
            }
        }
    }

    if(prms$mate.system == 'noSS.MPM'){
        if(male.biased){
            Nf.in <- length(f.ind)
            mated.f <- rbinom(Nf.in, 1, prob=(ff[f.ind,'pheno.sexf']/prms$L))
            mated.f.inds <- f.ind[mated.f==1]
            names(mated.f.inds) <- NULL
            n.matings <- length(mated.f.inds)
            mated <- vector(length = n.matings)

            avail.mates <- m.ind
            if(length(avail.mates)==1){
              mated[] <- avail.mates  
            }else{
            mated[] <- sample(avail.mates, n.matings,
                            replace=FALSE)
            }
            moms <- mated.f.inds
            dads <- mated
        }else{
            Nm.in <- length(m.ind)
            mated.m <- rbinom(Nm.in, 1, prob=(mm[m.ind,'pheno.sexm']/prms$L))
            mated.m.inds <- m.ind[mated.m==1]
            names(mated.m.inds) <- NULL
            n.matings <- length(mated.m.inds)
            mated <- vector(length = n.matings)

            avail.mates <- f.ind
            if(length(avail.mates)==1){
              mated[] <- avail.mates  
            }else{
            mated[] <- sample(avail.mates, n.matings,
                            replace=FALSE)
            }
            moms <- mated
            dads <- mated.m.inds
        }  
    }

    if(prms$mate.system == 'CTM'){
        n.matings <- pmin(length(f.ind), length(m.ind))
        mated <- vector(length = n.matings)
        if(length(f.ind)==1 & length(m.ind)==1){
            moms <- f.ind
            dads <- m.ind
        }else{
            if(male.biased){
                avail.mates <- m.ind
                mated[] <- sample(avail.mates, n.matings,
                                    prob = (mm[avail.mates,'pheno.sexm']+0.0001),
                                    replace=FALSE)
                moms <- f.ind
                dads <- mated
            }else{
                avail.mates <- f.ind
                mated[] <- sample(avail.mates, n.matings,
                                    prob = (ff[avail.mates,'pheno.sexf']+0.0001),
                                    replace=FALSE) 
                moms <- mated
                dads <- m.ind
            }
        }    
    }
  
    if(prms$mate.system == 'MPM'){
        if(male.biased){
            Nf.in <- length(f.ind)
            mated.f <- rbinom(Nf.in, 1, prob=(ff[f.ind,'pheno.sexf']/prms$L))
            mated.f.inds <- f.ind[mated.f==1]
            names(mated.f.inds) <- NULL
            n.matings <- length(mated.f.inds)
            mated <- vector(length = n.matings)

            avail.mates <- m.ind
            if(length(avail.mates)==1){
              mated[] <- avail.mates  
            }else{
            mated[] <- sample(avail.mates, n.matings,
                            prob = (mm[avail.mates,'pheno.sexm']+0.0001),
                            replace=FALSE)
            }
            moms <- mated.f.inds
            dads <- mated
        }else{
            Nm.in <- length(m.ind)
            mated.m <- rbinom(Nm.in, 1, prob=(mm[m.ind,'pheno.sexm']/prms$L))
            mated.m.inds <- m.ind[mated.m==1]
            names(mated.m.inds) <- NULL
            n.matings <- length(mated.m.inds)
            mated <- vector(length = n.matings)

            avail.mates <- f.ind
            if(length(avail.mates)==1){
              mated[] <- avail.mates  
            }else{
            mated[] <- sample(avail.mates, n.matings,
                            prob = (ff[avail.mates,'pheno.sexf']+0.0001),
                            replace=FALSE)
            }
            moms <- mated
            dads <- mated.m.inds
        }  
    }

## parent info
  parent.ids <- cbind(mom.id = ff[moms,'ID'],
                      dad.id = mm[dads,'ID'],
                      mom.ind = moms,
                      dad.ind = dads) 
  return(parent.ids)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Mutation ~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to mutate the loci at probability 'mu'. Since each loci is
## just a 0 or 1, I just need to determine if they switched and swap
## the ones that do. 

mutation <- function(loci, prms){
    if(prms$mu == 0){return(loci)}
  mut <- rbinom(length(loci), 1, prob=prms$mu)
  loci[mut>0] <- loci[mut>0] == FALSE

  return(loci)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~ Maturation ~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to update the adult population when there are offspring
## that have matured. When a juvenile matures, I assign it its traits
## through random, independent assorting of the parental genotypes. We
## also assign a sex here. 

update.adults <- function(pop=pop, prms=prms, coffin=coffin, matured, gen){
    n.matured <- nrow(matured)
    new.adults <- array(NA, dim = c(n.matured, ncol(pop$males)))
    dimnames(new.adults)[[2]] <- dimnames(pop$males)[[2]]
    
 
    p.traits <- array(rbinom(prms$L * prms$num.traits * n.matured, 1, prob=0.5) + 1,
                      dim=c(n.matured, (prms$L * prms$num.traits)))

    rand.assort <- array(rbinom(prms$L * prms$num.traits * n.matured, 1, prob=0.5) + 1,
                      dim=c(n.matured, (prms$L * prms$num.traits)))
    p.traits <- rand.assort                  
 
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Random, independent assortment ~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    new.adults[] <- t(sapply(1:nrow(matured), function(ii){
                par.traits <- array(NA, dim=c((prms$L * prms$num.traits), 2))
                
                if(matured[ii,'mom.id'] %in% pop[['females']][,'ID']){
                    par.traits[,1] <- pop[['females']][pop[['females']][,'ID'] %in% matured[ii,'mom.id'],1:(prms$L * prms$num.traits)]
                }else{
                    par.traits[,1] <- coffin[['dead.f']][coffin[['dead.f']][,'ID'] %in% matured[ii,'mom.id'],1:(prms$L * prms$num.traits)]  
                }

                if(matured[ii,'dad.id'] %in% pop[['males']][,'ID']){
                    par.traits[,2] <- pop[['males']][pop[['males']][,'ID'] %in% matured[ii,'dad.id'],1:(prms$L * prms$num.traits)]
                }else{
                    par.traits[,2] <- coffin[['dead.m']][coffin[['dead.m']][,'ID'] %in% matured[ii,'dad.id'],1:(prms$L * prms$num.traits)] 
                }
                traits <- sapply(1:(prms$L * prms$num.traits), function(xx) par.traits[xx,p.traits[ii,xx]])   
                traits <- c(traits, rbinom(1,1,prms$birth.ratio), 0, 1, matured[ii,'patch'], 0, 0, 0, 0)

                traits
            }
        )
    )

  ## ~~~~~~~~~~~~~~~
  ## ~~ Mutate ~~~~~
  ## ~~~~~~~~~~~~~~~
    loci <- 1:(prms$L*4)
    if(prms$ss.off == T){
        trait.set <- (loci-(loci%%prms$L))/prms$L
        trait.inds <- which(trait.set == 1 | trait.set == 3) + 1
        
        new.adults[,trait.inds] <- mutation(loci = new.adults[,trait.inds], prms)
    }else{
        new.adults[,loci] <- mutation(loci = new.adults[,loci], prms)
    }

    new.females <- new.adults[which(new.adults[,'sex']==0),,drop=FALSE]
    new.males <- new.adults[which(new.adults[,'sex']==1),,drop=FALSE]

    natal.patches.f <- new.females[,'patch'] 
    natal.patches.m <- new.males[,'patch']

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Offspring dispersal ~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(prms$n.patches == 2){
    for(disp in 1:prms$n.patches){
     f.disp.inds <- rbinom(nrow(new.females), 1, prob=prms$f.disp)==1
     m.disp.inds <- rbinom(nrow(new.males), 1,prob=prms$m.disp)==1
     
     new.females[f.disp.inds,'patch'] <- (new.females[f.disp.inds,'patch'] %%2)+1
     new.males[m.disp.inds,'patch'] <- (new.males[m.disp.inds,'patch'] %%2)+1
   } 
  } 
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Assign Phenotypes ~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~

    trait.set <- (loci-(loci%%prms$L))/prms$L

    new.females[,'pheno.sexm'] <- rowSums(new.females[,c(1,which(trait.set==0)+1), drop=FALSE])
    new.females[,'pheno.pcm']  <- rowSums(new.females[,which(trait.set==1)+1, drop=FALSE])
    new.females[,'pheno.sexf'] <- rowSums(new.females[,which(trait.set==2)+1, drop=FALSE])
    new.females[,'pheno.pcf']  <- rowSums(new.females[,which(trait.set==3)+1, drop=FALSE])
    
    new.males[,'pheno.sexm'] <- rowSums(new.males[,c(1,which(trait.set==0)+1), drop=FALSE])
    new.males[,'pheno.pcm']  <- rowSums(new.males[,which(trait.set==1)+1, drop=FALSE])
    new.males[,'pheno.sexf'] <- rowSums(new.males[,which(trait.set==2)+1, drop=FALSE])
    new.males[,'pheno.pcf']  <- rowSums(new.males[,which(trait.set==3)+1, drop=FALSE])

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Assign Individual ID ~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  new.females[,'ID'] <- (pop$id.counter['female.id'] + 1):
    (pop$id.counter['female.id'] + nrow(new.females)) 
  pop$id.counter['female.id'] <- pop$id.counter['female.id'] +
    nrow(new.females)
    
  new.males[,'ID'] <- (pop$id.counter['male.id'] + 1):
                      (pop$id.counter['male.id'] + nrow(new.males)) 
  pop$id.counter['male.id'] <- pop$id.counter['male.id'] + nrow(new.males) 

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Append to population ~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  pop$males <- rbind(pop$males, new.males)
  pop$females <- rbind(pop$females, new.females)

  if(length(new.males > 0)){
    new.males <- cbind(new.males, matings=0, off.surv=0,
                       patch.natal=natal.patches.m, total.off=0, lifespan=0)
  }else{
    new.males <- cbind(new.males, matings=integer(0), off.surv=integer(0),
                       patch.natal=integer(0), total.off=integer(0), lifespan=integer(0))
  }

  if(length(new.females > 0)){
    new.females <- cbind(new.females, matings=0, off.surv=0,
                         patch.natal=natal.patches.f, total.off=0, lifespan=0)
  }else{
    new.females <- cbind(new.females, matings=integer(0), off.surv=integer(0),
                       patch.natal=integer(0), total.off=integer(0), lifespan=integer(0))
  }
  
  updates <- list(pop=pop,
                  new.males=new.males[,c('sex', 'ID',
                                         'patch',
                                         'pheno.pcm',
                                         'pheno.sexm',
                                         'matings', 'off.surv',
                                         'patch.natal', 'total.off','lifespan')],
                  new.females=new.females[,c('sex', 'ID',
                                             'patch',
                                             'pheno.pcf',
                                             'pheno.sexf',
                                             'matings', 'off.surv',
                                             'patch.natal','total.off','lifespan')]) 
  return(updates)
}


update.fitbit <- function(fitbit, parents, broods, gen=gen){
    
    idm.inds <- match(parents[,'dad.id'], fitbit$fit.m[,'ID'])
    idf.inds <- match(parents[,'mom.id'], fitbit$fit.f[,'ID'])
    m.inds <- !is.na(idm.inds)
    f.inds <- !is.na(idf.inds)

    fitbit$fit.m[idm.inds[m.inds],'total.off'] <- fitbit$fit.m[idm.inds[m.inds],'total.off'] + broods[m.inds] 
    fitbit$fit.f[idf.inds[f.inds] ,'total.off']<- fitbit$fit.f[idf.inds[f.inds],'total.off'] + broods[f.inds] 
    fitbit$fit.f[idf.inds[f.inds],'matings'] <-  fitbit$fit.f[idf.inds[f.inds],'matings'] + 1
    fitbit$fit.m[idm.inds[m.inds],'matings'] <-  fitbit$fit.m[idm.inds[m.inds],'matings'] + 1
    return(fitbit)
}

update.coffin <- function(surv, foc.pop, foc.dead, juvs, foc.pop.id){
    if(length(surv) > 0){
        foc.dead <- rbind(foc.dead[,,drop=FALSE], foc.pop[surv == 0,,drop=FALSE])
    }
    foc.dead <- foc.dead[foc.dead[,'ID'] %in% juvs[,foc.pop.id],,drop=FALSE] 
    return(foc.dead)
}

update.lifespan <- function(pop, fitbit){
    m.inds <- match(pop$males[,'ID'], fitbit$fit.m[,'ID'])
    f.inds <- match(pop$females[,'ID'], fitbit$fit.f[,'ID']) 
    if(length(m.inds)>0){
        fitbit$fit.m[m.inds,'lifespan'] <- fitbit$fit.m[m.inds,'lifespan'] + 1
    }
    if(length(f.inds)>0){
        fitbit$fit.f[f.inds,'lifespan'] <- fitbit$fit.f[f.inds,'lifespan'] + 1
    }
    
    return(fitbit)
}