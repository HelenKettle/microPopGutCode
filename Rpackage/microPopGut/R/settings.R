#Model settings

settings=list(
    
    nbox=3, #number of model compartments 
    time.step=1/24/60, #n*one minute

    numDays=0.01,#,
    seed=1, #used for stochastic meal comp only (fluc.subst.comp=TRUE)
    fix.pH=TRUE,
    fluc.inflow=TRUE,
    fluc.subst.comp=TRUE,
    plotInflow=FALSE,
    showPlots=TRUE,
    BMpd=0, #number of bowel movements per day
    numStrains=1,
    startFromFile=1,#0,
    startFileName='numDays7seed1fix.pHFALSEfluc.inflowTRUEfluc.subst.compTRUEnumStrains1BMpd0C.from.diet50P.from.diet10W.from.diet1100transitTime1.25runNumber1RS.frac0.68LPc7gamma.mag1',
    previousNumDays=NA,
    RS.frac=0.78,#0.78, #i.e. RS/(RS+NSP) assuming RS+NSP=C
    RS.frac.mucin=0.5,
    
    #meal times (24h clock)
    breakfast.start=7,
    lunch.start=13,
    dinner.start=19,
    meal.duration.h=0.5, #(in hours)
    time.to.reach.colon.h=7, #time through stomach and small intestine etc (h)

    load.microbe.files=FALSE,
    microbesFolder='Groups/',
    save.pH=TRUE,
    
    pHLimit=TRUE,
    
    SCFAnames=c('Acetate','Butyrate','Propionate'),
    carbs=c('RS','NSP'),
    waterName='H2O', #'W'

    start.BM.time=list(7,c(7,19),c(7,15,21)), #time of bowel movement for 1, 2, or 3 BM per day
    
    microbeNames = c('Bacteroides','ButyrateProducers1','ButyrateProducers2','ButyrateProducers3','LactateProducers','PropionateProducers','Methanogens','NoButyFibreDeg','NoButyStarchDeg','Acetogens'),
    
    mucin.degrader=c('Bacteroides','NoButyFibreDeg'),

    relate.time.to.reach.colon.to.Tt=FALSE,

    modifyLactateProducers=TRUE
)
