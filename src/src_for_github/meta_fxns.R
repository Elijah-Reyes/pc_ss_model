

run.model <- function(...){
    prms <- make.prms(...)
    pop <- make.pop(prms)
    juveniles <- make.juv()
    coffin <- make.coffin(pp = pop$males)
    fitbit <- make.fitbit(pop=pop, prms=prms)
    pop.initial <- pop
    trait.sets <- cbind(c('Parental care', 'Sexual trait', 'Neutral trait'), 
                        c('care', 'sex', 'neutral'))
    
    if(prms$store.growth==TRUE){
        pop.growth <- array(NA, dim=c(4,prms$num.gens), dimnames=list(c('births', 'offspring.matured', 'offspring.deaths', 'adult.deaths')))
    }else{
        pop.growth <- 'Try running with store.growth==TRUE'
    }

    if(prms$store.hist==TRUE){
        hist.time <- array(NA, dim=c(11, prms$num.gens, 5, 3), 
                        dimnames=list(bins=seq(0,10,length.out=11), gen=1:prms$num.gens, 
                                        trait=c('pheno.pcm', 'pheno.pcf', 'pheno.sexm', 'pheno.sexf','neutral'), 
                                        patch=c('1','2','all')))
        brks <- seq(0,10,length.out=11)
    
        mk.dense <- function(xx) xx/sum(xx)
    }else{
        hist.time <- 'Try running with store.hist==TRUE'
    }

    out <- list(pop=pop, juveniles=juveniles, coffin=coffin, fitbit=fitbit, pop.growth=pop.growth)
    extinct.status <- 0
    pop.stats <- array(NA, dim=c(18, prms$num.gens, 3),
                       dimnames=list(stats = c('pop.size', 'osr',
                                               'male.pop', 'female.pop',
                                               'male.in.pop','female.in.pop',
                                               'male.care', 'female.care',
                                               'male.care.var', 'female.care.var', 
                                               'male.sex', 'female.sex',
                                               'male.sex.var', 'female.sex.var',
                                               'male.neutral', 'female.neutral',
                                               'male.neutral.var', 'female.neutral.var'),
                                     gen = 1:prms$num.gens,
                                     patch = c(1,2,'all')))

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~ Generation loop ~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(gen in 1:prms$num.gens){
    if(gen == 1){
        out$fitbit <- make.fitbit(pop=out$pop, prms=prms)
    }
    
    out <- lifecycle(pop=out$pop,
                     juveniles=out$juveniles,
                     coffin=out$coffin,
                     prms=prms,
                     fitbit=out$fitbit,
                     pop.growth=out$pop.growth,
                     gen=gen)
   
    for(pp in c(1:prms$n.patches, 'all')){
        mm <- out$pop$males[,,drop=FALSE]
        ff <- out$pop$females[,,drop=FALSE]
        if(pp != 'all'){
            foc.m <- mm[mm[,'patch']==pp, ]
            foc.f <- ff[ff[,'patch']==pp, ]
        }else{
            foc.m <- mm
            foc.f <- ff
        }

        pop.stats['pop.size',gen,pp] <- sum(nrow(foc.f), nrow(foc.m))
        if(pop.stats['pop.size',gen,pp] < 5){break}
        pop.stats['male.in.pop',gen,pp] <- sum(foc.m[,'state']==1)
        pop.stats['female.in.pop',gen,pp] <- sum(foc.f[,'state']==1)
        pop.stats['osr',gen,pp] <- pop.stats['male.in.pop',gen,pp]/pop.stats['pop.size',gen,pp]
        pop.stats['male.pop',gen,pp] <- nrow(foc.m)
        pop.stats['female.pop',gen,pp] <- nrow(foc.f)
        pop.stats['male.care',gen,pp] <- mean(foc.m[,'pheno.pcm'])
        pop.stats['female.care',gen,pp] <- mean(foc.f[,'pheno.pcf'])
        pop.stats['male.care.var',gen,pp] <- var(foc.m[,'pheno.pcm'])
        pop.stats['female.care.var',gen,pp] <- var(foc.f[,'pheno.pcf'])
        pop.stats['male.sex',gen,pp] <- mean(foc.m[,'pheno.sexm'])
        pop.stats['female.sex',gen,pp] <- mean(foc.f[,'pheno.sexf'])
        pop.stats['male.sex.var',gen,pp] <- var(foc.m[,'pheno.sexm'])
        pop.stats['female.sex.var',gen,pp] <- var(foc.f[,'pheno.sexf'])
        pop.stats['male.neutral',gen,pp] <-  mean(rowSums(foc.m[,sprintf('locus.%s', 41:50),drop=FALSE])) 
        pop.stats['female.neutral',gen,pp] <- mean(rowSums(foc.f[,sprintf('locus.%s', 41:50),drop=FALSE]))
        pop.stats['male.neutral.var',gen,pp] <-  var(rowSums(foc.m[,sprintf('locus.%s', 41:50),drop=FALSE])) 
        pop.stats['female.neutral.var',gen,pp] <- var(rowSums(foc.f[,sprintf('locus.%s', 41:50),drop=FALSE]))

        if(prms$store.hist==TRUE){
            hist.time[,gen,'pheno.pcm',pp] <- mk.dense(table(factor(foc.m[,'pheno.pcm'], level=0:10)))
            hist.time[,gen,'pheno.pcf',pp] <- mk.dense(table(factor(foc.f[,'pheno.pcf'], level=0:10))) 
            hist.time[,gen,'pheno.sexm',pp] <- mk.dense(table(factor(foc.m[,'pheno.sexm'], level=0:10))) 
            hist.time[,gen,'pheno.sexf',pp] <- mk.dense(table(factor(foc.f[,'pheno.sexf'], level=0:10))) 
            hist.time[,gen,'neutral',pp] <- mk.dense(table(factor(c(rowSums(foc.m[,sprintf('locus.%s', 41:50), drop=FALSE]),
                                                        rowSums(foc.f[,sprintf('locus.%s', 41:50), drop=FALSE])), 
                                                    level=0:10))) 
        }
    }
    if(1 %in% (pop.stats['pop.size',gen,] < 30)){
        extinct.status <- 1
        break
    }
    
## Plots results as you run to watch the pop over time    
     if(prms$run.watch == T){
      if(gen%%prms$print.every == 0){
        cat(sprintf('gen:%s  pop.size.1:%s pop.size.2:%s osr(m:f):%2.3f f\\mcare:%2.3f\\%2.3f f\\msex:%2.3f\\%2.3f\n',      
                    gen, pop.stats['pop.size', gen, 1],
                    pop.stats['pop.size', gen, 2],
                    pop.stats['osr', gen, 1],
                    pop.stats['female.care', gen, 1],
                    pop.stats['male.care', gen, 1],
                    pop.stats['female.sex', gen, 1],
                    pop.stats['male.sex', gen, 1]))}
        if(gen%%500 == 0){
          m <- matrix(1:6, ncol=3)  
          layout(m)
          par(mar=c(1,4,0.5,1), oma=c(3,1,2,1))

          plot.run(y.lim=c(0,(1.25 * max(pop.stats['pop.size', ,1], na.rm=T))), 
                   stat='Population size', num.gens=prms$num.gens)

          for(patch in 1:prms$n.patches){
            lines(x=1:prms$num.gens, y=pop.stats['pop.size',,patch],
                  type='l', col='black', lty=patch)
            add.lines(n.gens=prms$num.gens, stat=pop.stats, trait='pop', patch)      
          }
          axis(side=2, at=seq(0, 1000, length.out=6),
               labels=seq(0, 1000, length.out=6))
          
          for(gg in 1:3){
                plot.run(y.lim=c(0,10), stat=trait.sets[gg,1], num.gens=prms$num.gens)
                for(patch in 1:prms$n.patches){
                    add.lines(n.gens=prms$num.gens, stat=pop.stats, 
                              trait=trait.sets[gg,2], patch)
                }
                trait.axis()
          }
          axis(side=1, at=seq(0,prms$num.gens,length.out=6),
               labels=seq(0,prms$num.gens,length.out=6), tck=0.01)
          mtext('Generation', side=1, at=prms$num.gens/2, line=2)     

        }
     }
  }
  
  out$fitbit$fit.m <- out$fitbit$fit.m[which(out$fitbit$fit.m[,'patch.natal']>0),]
  out$fitbit$fit.f <- out$fitbit$fit.f[which(out$fit$fit.f[,'patch.natal']>0),]

  out$fitbit$fit.stats['male',1,1,] <- colMeans(out$fitbit$fit.m[which(out$fitbit$fit.m[,'patch']==1 & 
                                                                   out$fitbit$fit.m[,'patch.natal']==1),
                                                                   c('matings', 'off.surv', 'total.off')]) 
  out$fitbit$fit.stats['male',1,2,] <- colMeans(out$fitbit$fit.m[which(out$fitbit$fit.m[,'patch']==1 & 
                                                                   out$fitbit$fit.m[,'patch.natal']==2),
                                                                   c('matings', 'off.surv', 'total.off')])  
  out$fitbit$fit.stats['male',2,1,] <- colMeans(out$fitbit$fit.m[which(out$fitbit$fit.m[,'patch']==2 & 
                                                                   out$fitbit$fit.m[,'patch.natal']==1),
                                                                   c('matings', 'off.surv', 'total.off')])    
  out$fitbit$fit.stats['male',2,2,] <- colMeans(out$fitbit$fit.m[which(out$fitbit$fit.m[,'patch']==2 & 
                                                                   out$fitbit$fit.m[,'patch.natal']==2),
                                                                   c('matings', 'off.surv', 'total.off')])   
  out$fitbit$fit.stats['female',1,1,] <- colMeans(out$fitbit$fit.f[which(out$fitbit$fit.f[,'patch']==1 & 
                                                                   out$fitbit$fit.f[,'patch.natal']==1),
                                                                   c('matings', 'off.surv', 'total.off')]) 
  out$fitbit$fit.stats['female',1,2,] <- colMeans(out$fitbit$fit.f[which(out$fitbit$fit.f[,'patch']==1 & 
                                                                   out$fitbit$fit.f[,'patch.natal']==2),
                                                                   c('matings', 'off.surv', 'total.off')])  
  out$fitbit$fit.stats['female',2,1,] <- colMeans(out$fitbit$fit.f[which(out$fitbit$fit.f[,'patch']==2 & 
                                                                   out$fitbit$fit.f[,'patch.natal']==1),
                                                                   c('matings', 'off.surv', 'total.off')])    
  out$fitbit$fit.stats['female',2,2,] <- colMeans(out$fitbit$fit.f[which(out$fitbit$fit.f[,'patch']==2 & 
                                                                   out$fitbit$fit.f[,'patch.natal']==2),
                                                                   c('matings', 'off.surv', 'total.off')])      
  pop <- out$pop                                                                                                                                                                                                                                                                                                                            
  fitness <- out$fitbit
  growth <- out$pop.growth
  return(list(final.pop=pop, prms=prms, pop.stats=pop.stats, 
              extinct.status=extinct.status, pop.initial=pop.initial, 
              fitness=fitness, pop.growth=growth, hist.time=hist.time))
}


plot.run <- function(y.lim, stat, num.gens=prms$num.gens){
    plot(x = NA, y= NA, yaxt='n', xaxt='n', 
         xlim=c(1, num.gens), ylim=y.lim,
         xlab='', ylab=stat)
}

add.lines <- function(n.gens=prms$num.gens, stat=pop.stats, trait, patch){
    ifelse(patch==1, clrs <- c('red', 'blue'), clrs <- c('coral1','cyan3'))
    lines(x=1:n.gens, y=stat[sprintf('female.%s', trait),,patch],
                    type='l', col=clrs[1], lty=patch)
    lines(x=1:n.gens, y=stat[sprintf('male.%s', trait),,patch],
                    type='l', col=clrs[2], lty=patch)
}

trait.axis <- function(){
    axis(side=2, at=seq(0, 10, length.out=6),
         labels=seq(0, 10, length.out=6), tck=0.01)
}

## Checks to see which files did not run after finishing
output.checker <- function(index, save.dir){
    inds <-  1:nrow(index)
    saved.files <- list.files(save.dir)
    saved.rows <- gsub('.RData','',saved.files)

    run.files <- which(as.character(inds) %in% saved.rows)
    unrun.files <- which((as.character(inds) %in% saved.rows) == FALSE)

    message(length(run.files),'/',nrow(index), ' expected runs completed.\n',
            'To see which files failed, type "failed.runs"')
    failed.runs <- unrun.files   
    return(failed.runs)     
}

pkg.save <- function(data, save.every){
    gens.to.store <- seq(0,100000,by=save.every)
    gens.to.store[1] <- 1

    data$hist.time <- data$hist.time[,gens.to.store,,]
    data$pop.stats <- data$pop.stats[,gens.to.store,]

    data <- list(data$extinct.status, data$fitness, 
                 data$hist.time, data$pop.growth,  
                 data$pop.stats, data$prms) 
    return(data)
}