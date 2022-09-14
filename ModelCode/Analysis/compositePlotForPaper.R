rm(list=ls())
graphics.off()

#composite plot barplots fot Tt, BMpd, fix.pH for PNAS paper
#makes 3 plots which need to be put side by side


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

exptNamesVec=c('transitTime','BMpd','fix.pH')

for (lv in c(2,3,1)){
    
lead.var=exptNamesVec[lv]

numDays=28
box=3

bac.cols=brewer.pal(10,'Paired')
scfa.cols=brewer.pal(3,'Set1')
SCFAnames=c('Acetate','Butyrate','Propionate')
waterName='H2O'

batch.par.names=c('numDays','seed','fix.pH','fluc.inflow','fluc.subst.comp','numStrains','BMpd','C.from.diet','P.from.diet','W.from.diet','transitTime','runNumber','startFromFile','previousNumDays','RS.frac','LPc','gamma.mag')

microbeNames.short=c('B','BP1','BP2','BP3','LP','PP','M','NBFD','NBSD','A')

#max for plot lims
#if (lead.var=='transitTime'){
    scfa.ymax=140
    bac.ymax=40
#}else{
#    scfa.ymax=130
#    bac.ymax=15
#}

labels=c('RS.frac'='RS fraction','P.from.diet'='Dietary protein (g/d)','transitTime'='Transit Time (d)',
    'BMpd'='#bowel movements/d','fix.pH'='fixed pH','fluc.inflow'='Meals',C.from.diet='C from diet (g/d)')




exptFolder=exptNamesVec
names(exptFolder)=exptNamesVec

folderName=paste0(resFolder,exptFolder[lead.var],'/')

all.files=list.files(folderName)
Lf=length(all.files)

print(paste0('there are ',Lf,' files'))

ft=makeFileTable(all.files,batch.par.names)

ii=ft[,'numDays']==numDays


#find the files I want
fix.var=c('numDays'=paste(numDays)
          )

keep=NULL
for (i in 1:Lf){
    if (identical(ft[i,names(fix.var)],fix.var)){
        keep=c(keep,i)
    }
}

        
num.seeds=max(as.numeric(ft[keep,'seed']))

print(ft[keep,])

good.files.index=keep
Tb1=ft[keep,]
good.filenames=all.files[keep]
Lfgood=length(good.filenames)

num.lead.var.vals=length(unique(Tb1[,lead.var]))
lead.var.vals=unique(Tb1[,lead.var])

num.other.var.vals=Lfgood/num.lead.var.vals

bac.mat=array(NA,dim=c(Lfgood,10))
scfa.mat=array(NA,dim=c(Lfgood,3))
pH=array(NA,dim=c(Lfgood,3))
TSCFA=array(NA,dim=c(Lfgood,3))


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
        main.txt=paste0(labels[lead.var],' = ',inf[lead.var],', meals (seed=',inf['seed'],')')
    }else{
        main.txt=paste0(labels[lead.var],' = ',inf[lead.var],', continuous inflow')
    }
    
    
    bac.mat[f,]=colMeans(out$microbe.array[ii,,box],na.rm=TRUE)
    scfa.mat[f,]=colMeans(out$SCFA.array[ii,SCFAnames,box],na.rm=TRUE)
    for (b in 1:3){
        pH[f,b]=mean(out$ph.mat[ii,b],na.rm=TRUE)
        TSCFA[f,b]=sum(colMeans(out$SCFA.array[ii,SCFAnames,b]))
    }
    
}




#Plots of distal for bac and scfa
#average over seeds------------------------------------------
w=2
dev.new(width=w*2,height=2*3.5)

pct=1

par(mfrow=c(2,w),mar=c(5,4,3,0))


for (p in 1:2){
        
    for (fluc in c('TRUE','FALSE')){            

        if (p==1){
            nc=10
            mat=bac.mat
            cols=bac.cols
            ytxt='Biomass (g)'
            legn=microbeNames.short
            ymax=bac.ymax
            legsize=0.8
            legpos='topleft'
        }else{
            nc=3
            mat=scfa.mat
            cols=scfa.cols
            ytxt='SCFA (mM)'
            legn=SCFAnames
            ymax=scfa.ymax
            legsize=0.9
            legpos='topright'
        }
        
        if (fluc){
            main.txt=paste('meals')
        }else{
            main.txt='cts inflow'
        }
        avmat=matrix(NA,nrow=length(lead.var.vals),ncol=nc)
        rownames(avmat)=paste(lead.var.vals)
        
        for (ct in 1:length(lead.var.vals)){
            
            ii=(Tb1[,lead.var]==lead.var.vals[ct] & Tb1[,'fluc.inflow']==fluc)
            
            if (sum(ii)>1){
                avmat[ct,]=colMeans(mat[ii,])
            }else{
                avmat[ct,]=mat[ii,]
            }
        }
        
            #sort into numerical order
        if (nrow(avmat)>2){
            iii=sort(as.numeric(rownames(avmat)),index.return=TRUE)$ix
        }else{
            iii=seq(1,nrow(avmat))
        }
        
        
        if (lead.var=='transitTime'){
                    xlabel='Transit time (d)'
                    names=rownames(avmat[iii,])
        }else if (lead.var=='BMpd'){
                    names=c('cts','2 BM/d')
                    xlabel='Outflow'
        }else if (lead.var=='fix.pH'){
                    names=c('Vary','Fixed')
                    xlabel='pH'
        }    
             barplot(t(avmat[iii,]),col=cols,
                     main=main.txt,
                     names.arg=names,
                     ylab=ytxt,
                     xlab=xlabel,
                     ylim=c(0,ymax),
                     cex.lab=1.1,
                     cex.axis=1.1,
                     cex.names=1.1)
        legend(legpos,col=cols,legend=legn,pch=15,cex=legsize,bty='n')

        pct=pct+1

    #   if (lead.var=='transitTime'){
    #    if (pct%%3==0){
    #        blankPlot()
    #        legend('left',col=cols,legend=legn,pch=15,cex=1.6)
    #        pct=pct+1
    #    }
    #}
            
        
    }
}

if (lead.var=='fix.pH'){lead.var='fixpH'}
dev.copy2eps(file=paste0(lead.var,'barSmall.eps'))

}
