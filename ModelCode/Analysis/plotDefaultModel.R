rm(list=ls())
graphics.off()

#using base R graphics

library('RColorBrewer')
library('xtable')
library('abind')


source('makeFileTable.R')
source('loadOutput.R')
source('computeSCFA.R')
source('blankPlot.R')

resFolder='../../ModelRuns/'


include.cts=TRUE


lead.var='transitTime'
numDays=28
box=3


plotTS=TRUE

bac.cols=brewer.pal(10,'Paired')
scfa.cols=brewer.pal(3,'Set1')
SCFAnames=c('Acetate','Butyrate','Propionate')
waterName='H2O'

batch.par.names=c('numDays','seed','fix.pH','fluc.inflow','fluc.subst.comp','numStrains','BMpd','C.from.diet','P.from.diet','W.from.diet','transitTime','runNumber','startFromFile','previousNumDays','RS.frac','LPc','gamma.mag')

microbeNames.short=c('B','BP1','BP2','BP3','LP','PP','M','NBFD','NBSD','A')

#max for plot lims
scfa.ymax=140
bac.ymax=30

labels=c('transitTime'='Transit Time (d)')


folderName=paste0(resFolder,lead.var,'/')

summary=list()

all.files=list.files(folderName)
Lf=length(all.files)

print(paste0('there are ',Lf,' files'))

ft=makeFileTable(all.files,batch.par.names)

#find the files I want
#only plot seeds 1 and 3

ii=ft[,'seed']%in%c('1','4') & ft[,'numDays']=='28' & ft[,'transitTime']%in%c('1','2')

Tb1=ft[ii,]
good.filenames=all.files[ii]
Lfgood=length(good.filenames)

num.lead.var.vals=length(unique(Tb1[,lead.var]))
lead.var.vals=unique(Tb1[,lead.var])

num.other.var.vals=Lfgood/num.lead.var.vals

bac.mat=array(NA,dim=c(Lfgood,10))
scfa.mat=array(NA,dim=c(Lfgood,3))
pH=array(NA,dim=c(Lfgood,3))
TSCFA=array(NA,dim=c(Lfgood,3))

#only plot every other file for Tt
np=2
if (plotTS){
    dev.new(height=np*3,width=(num.other.var.vals+1)*3)
    par(mfcol=c(np,num.other.var.vals+1),mar=c(5,5,3,0.5))
    dev.new(height=np*3,width=(num.other.var.vals+1)*3)
    par(mfcol=c(np,num.other.var.vals+1),mar=c(5,5,3,0.5))
}
for (f in 1:Lfgood){

    file=good.filenames[f]
    
    out=loadOutput(file,folderName,waterName,SCFAnames,substrateNames=c('Protein','RS','NSP'))

    time=out$time
    Lt=length(time)
    
    microbeNames=names(out$microbe.array[1,,1])

    #av over 7-numDays
    if (max(time>8)){
        ii=time>7
    }else{
        ii=time>=0
    }

    inf=Tb1[f,]
    if (inf['fluc.inflow']){
        main.txt=paste0('Transit time = ',inf[lead.var],' d, meals (seed=',inf['seed'],')')
    }else{
        main.txt=paste0('Transit time = ',inf[lead.var],' d, cts inflow')
    }
    
    if (plotTS){
 
        dev.set(2)
        matplot(time,out$microbe.array[,microbeNames,box],type='l',col=bac.cols,lwd=2,lty=1,
                ylab='Biomass (g)',xlab='time (d)',
                main=main.txt,
                cex.lab=1.5,cex.axis=1.4,cex.main=1.4
                )
        dev.set(3)
        matplot(time,out$SCFA.array[,SCFAnames,box],type='l',col=scfa.cols,lwd=2,lty=1,
                ylab='SCFA (mM)',xlab='time (d)',
                main=main.txt,
                cex.lab=1.5,cex.axis=1.4,cex.main=1.4
                )
        
    
 
    }

    
    bac.mat[f,]=colMeans(out$microbe.array[ii,,box],na.rm=TRUE)
    scfa.mat[f,]=colMeans(out$SCFA.array[ii,SCFAnames,box],na.rm=TRUE)
    for (b in 1:3){
        pH[f,b]=mean(out$ph.mat[ii,b],na.rm=TRUE)
        TSCFA[f,b]=sum(colMeans(out$SCFA.array[ii,SCFAnames,b]))
    }
    
}


dev.set(2)
for (i in 1:np){
    blankPlot()
    legend('left',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.5,bty='n')
}
dev.set(3)
for (i in 1:np){
    blankPlot()
    legend('left',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.5,bty='n')
}

dev.set(2)
dev.copy2eps(file=paste0(lead.var,'TS.eps'))
