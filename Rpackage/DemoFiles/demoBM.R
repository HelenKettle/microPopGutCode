#Bowel Movements
#constant inflow

rm(list=ls())
#graphics.off()


microbeGroups=c(1,2) #pick groups to simulate
sim.time.h=48 #time to simulate in hours


library(devtools)
pf='/home/helen/MicroPop/MicroPopGutPackage/microPopGut'
load_all(pf)

microbeNames = c('Bacteroides','ButyrateProducers1','ButyrateProducers2','ButyrateProducers3','LactateProducers','PropionateProducers','Methanogens','NoButyFibreDeg','NoButyStarchDeg','Acetogens')

microbeNames.short = c('Bacteroides'='B','ButyrateProducers1'='BP1','ButyrateProducers2'='BP2','ButyrateProducers3'='BP3','LactateProducers'='LP','PropionateProducers'='PP','Methanogens'='M','NoButyFibreDeg'='NBFD','NoButyStarchDeg'='NBSD','Acetogens'='A')


m.out=microPopGut(
    numDays=sim.time.h/24,
    time.step=1/24/60,
    transitTime=1.25,
    microbeNames=microbeNames[microbeGroups],
    microbeNames.short=microbeNames.short[microbeGroups],
    bowel.movements = list(BM.duration.h = 15/60,
        frac.distal.emptied = 0.9,
        BMpd = 3,
        start.BM.time = list(7, c(7, 19), c(7, 15, 21)))
)

time=m.out$solution[[1]][,'time']
verification(m.out,start.av=0.8*max(time),fin.av=max(time))
dev.new()
plotMPG(m.out)

dev.copy2eps(file=paste0('ResultDemoBM_',length(microbeGroups),'groups_time',sim.time.h/24,'.eps'))
