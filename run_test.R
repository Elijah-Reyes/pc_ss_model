## Set working directory here:
setwd()
setwd('~/Dropbox/Projects/R_fromhage_jennions_for_GitHub/')
## setwd('~/Dropbox/R_fromhage_jennions')

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


tester.arr <- array(NA, dim=c(2,6,5000), dimnames=list(patch=1:2,
                                                       trait =
                                                         c('male.care',
                                                           'female.care',
                                                           'male.sex',
                                                           'female.sex',
                                                           'male.neutral',
                                                           'female.neutral'),
                                                       gen=1:5000))

tester.arr[1,,] <- out$pop.stats[ c('male.care',
                                    'female.care',
                                    'male.sex',
                                    'female.sex',
                                    'male.neutral',
                                    'female.neutral'),,1]

tester.arr[2,,] <- out$pop.stats[ c('male.care',
                                    'female.care',
                                    'male.sex',
                                    'female.sex',
                                    'male.neutral',
                                    'female.neutral'),,2]


pdf('dispersalPostEvolution_parentalCare_plane.pdf', width=14, height=7)
m <- matrix(1:2, ncol=2)
layout(m)
par(cex=1.2, tck=0.01, mgp=c(3,0.5,0), oma=c(3,4,3,2), mar=c(1,0.75,0.75,1))
plot(x=NA, y=NA, xaxt='n', yaxt='n', xlab='', ylab='',
     xlim=c(0,10), ylim=c(0,10))
lines(x=tester.arr[1,'female.care',],
      y=tester.arr[1,'male.care',], col='red')
lines(x=tester.arr[2,'female.care',],
      y=tester.arr[2,'male.care',], col='blue')

abline(0, 1, lty=2)
axis(1, at=seq(0,10,2), lab=seq(0,10,2))
axis(2, at=seq(0,10,2), lab=seq(0,10,2))
mtext('Female parental care', side=1, at=5, line=2, cex=1.4)
mtext('Male parental care', side=2, at=5, line=2, cex=1.4)
mtext('Dispersal = 0.01', side=3, line=0.75, outer=T, cex=1.4)

plot(x=NA, y=NA, xaxt='n', yaxt='n', xlab='', ylab='',
     xlim=c(0,10), ylim=c(0,10))
lines(x=tester.arr[1,'female.neutral',],
      y=tester.arr[1,'male.neutral',], col='red')
lines(x=tester.arr[2,'female.neutral',],
      y=tester.arr[2,'male.neutral',], col='blue')

abline(0, 1, lty=2)
axis(1, at=seq(0,10,2), lab=seq(0,10,2))
## axis(2, at=seq(0,10,2), lab=seq(0,10,2))
mtext('Female neutral trait', side=1, at=5, line=2, cex=1.4)
mtext('Male neutral trait', side=2, at=5, line=0.5, cex=1.4)

dev.off()



xx <- out$final.pop$fitbit$fit.f
yy <- out$final.pop$fitbit$fit.f

p1.inds <- which(xx[,'patch']==1)
p2.inds <- which(xx[,'patch']==2)

plot(x=xx[p1.inds,'parCare'],y=yy[p1.inds,'off.surv'])









layout(matrix(1:4, ncol=2))

plot(x=out$fitness$fit.m[,'parCare'], y=out$fitness$fit.m[,'off.surv'])
abline(lm(out$fitness$fit.m[,'off.surv']~out$fitness$fit.m[,'parCare']), col='red')

plot(x=out$fitness$fit.f[,'parCare'], y=out$fitness$fit.f[,'off.surv'])
abline(lm(out$fitness$fit.f[,'off.surv']~out$fitness$fit.f[,'parCare']), col='red')

plot(x=foo$fitness$fit.m[,'parCare'], y=foo$fitness$fit.m[,'off.surv'])
abline(lm(foo$fitness$fit.m[,'off.surv']~foo$fitness$fit.m[,'parCare']), col='red')

plot(x=foo$fitness$fit.f[,'parCare'], y=foo$fitness$fit.f[,'off.surv'])
abline(lm(foo$fitness$fit.f[,'off.surv']~foo$fitness$fit.f[,'parCare']), col='red')
