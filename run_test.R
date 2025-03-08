setwd(##Your Project folder)
rm(list=ls())

source('src/init.R')
   
out <- run.model(male.care.inits = c(0,9),
                 female.care.inits = c(9,0),
                 male.sex.inits = c(7,1),
                 female.sex.inits = c(1,7),
                 m.disp=0.005,
                 f.disp=0.005,
                 n.patches = 2,
                 mate.system='noSS.CTM', ## See prms for possible choices
                 num.gens=2e4,
                 ss.off=FALSE,
                 store.fit=TRUE,
                 store.growth=TRUE,
                 store.hist=TRUE,
                 inits.fixed=TRUE,
                 run.watch=TRUE,
                 print.every=5000) 

