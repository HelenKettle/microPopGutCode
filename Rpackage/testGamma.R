rm(list=ls())
graphics.off()

#meals

#---------------------GUT PACKAGE -----------------------
library(devtools)
pf='/home/helen/MicroPop/MicroPopGutPackage/microPopGut'
load_all(pf)

num.days=7
breakfast.start=7
lunch.start=13
dinner.start=19
time.step.h=1/60
meal.duration.h=0.5
time.to.reach.colon.h=7
#mealComposition------------------------------
substrateNames=c('A','B','C')
substrateMeans=c(A=1,B=2,C=3)
scale=c(A=0.5,B=1,C=1.5)
#make substrate mat------------

RS.frac=0.78
makePlot=FALSE
waterName='water'
fluc.comp=TRUE
waterMean=1100

for (gamma.mag in c(0.1,1,1.9)){

    out=makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag,RS.frac,showPlot=TRUE,saveFig=TRUE,waterName,waterMean)

    system(paste0('mv mealComp.eps OfflineVignette/mealComp',gamma.mag,'.eps'))
    system(paste0('mv outflowFromSmallIntestine.eps OfflineVignette/outflowFromSmallIntestine',gamma.mag,'.eps'))
    
}

