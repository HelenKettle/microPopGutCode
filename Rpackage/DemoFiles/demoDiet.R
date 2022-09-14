#constant inflow

rm(list=ls())
#graphics.off()


microbeGroups=c(1,2) #pick groups to simulate
sim.time.h=24 #time to simulate in hours


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

    #initial mass in each compartment
    init=list(
        C=2, #carbohydrate (g)
        P=0, #protein (g)
        B=10, #biomass (g)
        Acetate=0.3606, #g
        Propionate=0.1482, #g
        Butyrate=0.1762, #g
        W=100), #water (g)

    #inflow from diet
    inflow=list(
        C=60, #carbohydrate (g/d)
        P=0, #protein (g/d)
        W=1100, #water (g/d)
        RS.frac=0.78) #fraction of C that is resistant starch (rest is NSP)
    
)

time=m.out$solution[[1]][,'time']
verification(m.out,start.av=0.8*max(time),fin.av=max(time))
dev.new()
plotMPG(m.out)

dev.copy2eps(file=paste0('ResultDemoDiet_',length(microbeGroups),'groups_time',sim.time.h/24,'.eps'))

