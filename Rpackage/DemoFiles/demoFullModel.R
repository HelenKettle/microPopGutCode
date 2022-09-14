rm(list=ls())
#graphics.off()

install.packages('../../microPopGut_1.0.tar.gz')
library(microPopGut)

#library(devtools)
#pf='/home/helen/MicroPop/MicroPopGutPackage/microPopGut'
#load_all(pf)


microbeGroups=1:10 #pick groups to simulate
sim.time.h=2 #time to simulate in hours



microbeNames = c('Bacteroides','ButyrateProducers1','ButyrateProducers2','ButyrateProducers3','LactateProducers','PropionateProducers','Methanogens','NoButyFibreDeg','NoButyStarchDeg','Acetogens')

microbeNames.short = c('Bacteroides'='B','ButyrateProducers1'='BP1','ButyrateProducers2'='BP2','ButyrateProducers3'='BP3','LactateProducers'='LP','PropionateProducers'='PP','Methanogens'='M','NoButyFibreDeg'='NBFD','NoButyStarchDeg'='NBSD','Acetogens'='A')

sim.time.h=28 #time to simulate in hours

m.out=microPopGut(
    numDays=sim.time.h/24,
    time.step=1/24/60,
    transitTime=1.25,
    microbeNames=microbeNames[microbeGroups],
    microbeNames.short=microbeNames.short[microbeGroups],
    meals=list(fluc.inflow=TRUE,seed=2,fluc.subst.comp=TRUE,plotInflow=TRUE,saveInflowFig=TRUE)
)

time=m.out$solution[[1]][,'time']
verification(m.out,start.av=0.8*max(time),fin.av=max(time))
plotMPG(m.out)

dev.copy2eps(file=paste0('ResultDemoFullModel_',length(microbeGroups),'groups_time',sim.time.h/24,'.eps'))
