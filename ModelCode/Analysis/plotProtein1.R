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

box=3

plotTS=TRUE

bac.cols=brewer.pal(10,'Paired')
scfa.cols=brewer.pal(3,'Set1')
SCFAnames=c('Acetate','Butyrate','Propionate')
waterName='H2O'

batch.par.names=c('numDays','seed','fix.pH','fluc.inflow','fluc.subst.comp','numStrains','BMpd','C.from.diet','P.from.diet','W.from.diet','transitTime','runNumber','startFromFile','previousNumDays','RS.frac','LPc','gamma.mag')

microbeNames.short=c('B','BP1','BP2','BP3','LP','PP','M','NBFD','NBSD','A')

lead.var='C.from.diet'

folderName=paste0(resFolder,'DuncanNew/')

all.files=list.files(folderName)
Lf=length(all.files)

print(paste0('there are ',Lf,' files'))

ft=makeFileTable(all.files,batch.par.names)

CfromdietVec=seq(10,50,10)

Tt='1.5'#paste(c(1,3))

RSfvec=c('0.2','0.78')

diets=paste(CfromdietVec)

bac.array.all=array(NA,dim=c(length(diets),10,2,length(RSfvec)),dimnames=list(diets,microbeNames.short,c('cts','meals'),paste(RSfvec)))
scfa.array.all=array(NA,dim=c(length(diets),3,2,length(RSfvec)),dimnames=list(diets,SCFAnames,c('cts','meals'),paste(RSfvec)))

numDays=28


for (RSf in RSfvec){


all.bac=array(NA,dim=c(length(diets),10,2),
        dimnames=list(diets,microbeNames.short,c('cts','av.seed')))

all.scfa=array(NA,dim=c(length(diets),3,2),
        dimnames=list(diets,SCFAnames,c('cts','av.seed')))

for (diet in diets){

#find the files I want
        fix.var=c('numDays'=paste(numDays),
            RS.frac=RSf,
            C.from.diet=diet
                  )

        keep=NULL
        for (i in 1:Lf){
            if (identical(ft[i,names(fix.var)],fix.var)){
                keep=c(keep,i)
            }
        }
    
        if (is.null(keep)){
            break('no files fit criteria')
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
        
                                        #only plot every other file for Tt
        np=num.lead.var.vals
        if (plotTS & diet==diets[1]){
            dev.new(height=length(diets)*2,width=(num.other.var.vals+1)*2)
            par(mfrow=c(length(diets),num.other.var.vals+1),c(4,5,2,0.5))
            dev.new(height=length(diets)*2,width=(num.other.var.vals+1)*2)
        par(mfrow=c(length(diets),num.other.var.vals+1),c(4,5,2,0.5))
        }
        for (f in 1:Lfgood){
            
            file=good.filenames[f]
            
            out=loadOutput(file,folderName,waterName,SCFAnames,substrateNames=c('Protein','RS','NSP'))
            
            time=out$time
            Lt=length(time)
            
            microbeNames=names(out$microbe.array[1,,1])

    #av over 7-numDays if more than 8 days, otherwise av over whole time
            if (max(time>8)){
                ii=time>7
            }else{
                ii=time>=0
            }
            
            inf=Tb1[f,]
            if (inf['fluc.inflow']){
                main.txt=paste0('C from diet = ',inf[lead.var],', meals (seed=',inf['seed'],')')
            }else{
                main.txt=paste0('C from diet= ',inf[lead.var],', continuous inflow')
            }
            
            if (plotTS){
                dev.set(2)
                matplot(time,out$microbe.array[,microbeNames,box],type='l',col=bac.cols,lwd=2,lty=1,
                        ylab='biomass (g)',xlab='time (d)',
                        main=main.txt,
                        cex.lab=1.2,cex.axis=1.2
                        )
                dev.set(3)
                matplot(time,out$SCFA.array[,SCFAnames,box],type='l',col=scfa.cols,lwd=2,lty=1,
                        ylab='SCFA (mM)',xlab='time (d)',
                        main=main.txt,
                        cex.lab=1.2,cex.axis=1.2
                        )
                
            }
            
            
            bac.mat[f,]=colMeans(out$microbe.array[ii,,box],na.rm=TRUE)
            scfa.mat[f,]=colMeans(out$SCFA.array[ii,SCFAnames,box],na.rm=TRUE)
            for (b in 1:3){
                pH[f,b]=mean(out$ph.mat[ii,b],na.rm=TRUE)
                TSCFA[f,b]=sum(colMeans(out$SCFA.array[ii,SCFAnames,b]))
            }
            
        }
        
        
        if (plotTS){
            dev.set(2)
            for (i in 1:np){
                blankPlot()
                legend('left',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.1)
            }
            dev.set(3)
            for (i in 1:np){
                blankPlot()
                legend('left',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.1)
            }
            
            dev.set(2)
            dev.copy2eps(file=paste0(lead.var,'TS.eps'))
        }
        
        ii.cts=Tb1[,'fluc.inflow']==FALSE
        ii.fluc=Tb1[,'fluc.inflow']==TRUE
        if (sum(ii.cts)>0){
            all.bac[diet,,'cts']=bac.mat[ii.cts,]
            all.scfa[diet,,'cts']=scfa.mat[ii.cts,]
        }
        if (sum(ii.fluc>0)){
            if (num.seeds>1){
                all.bac[diet,,'av.seed']=colMeans(bac.mat[ii.fluc,])
                all.scfa[diet,,'av.seed']=colMeans(scfa.mat[ii.fluc,])
            }else{
                all.bac[diet,,'av.seed']=bac.mat[ii.fluc,]
                all.scfa[diet,,'av.seed']=scfa.mat[ii.fluc,]
            }
        }
        
    } #diet

    bac.array.all[,,,paste(RSf)]=all.bac
    scfa.array.all[,,,paste(RSf)]=all.scfa

}#RSf 

#plot mean SCFA
dev.new(height=7,width=10)
par(mfrow=c(2,length(RSfvec)*2),mar=c(5,5,3,2))

for (RSf in RSfvec){
    barplot(t(bac.array.all[,,'cts',RSf]),
            col=bac.cols,
            names.arg=CfromdietVec,
            xlab='C inflow (g/d)',ylab='Biomass (g)',
            ylim=c(0,30),
             main=paste('RS frac ',RSf,'(cts inflow)'),
            cex.names=1.2,cex.axis=1.2,cex.lab=1.2,cex.main=1.4
                        )
    legend('topleft',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.2,bty='n')
    
    barplot(t(bac.array.all[,,'meals',RSf]),
            col=bac.cols,
            names.arg=CfromdietVec,
            xlab='C inflow (g/d)',ylab='Biomass (g)',
            ylim=c(0,30),
            main=paste('RS frac ',RSf,' (meals)'),
            cex.names=1.2,cex.axis=1.2,cex.lab=1.2,cex.main=1.4
            )
    legend('topleft',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.2,bty='n')
    

}

for (RSf in RSfvec){
    barplot(t(scfa.array.all[,,'cts',RSf]),
            names.arg=CfromdietVec,
            col=scfa.cols,
            xlab='C inflow (g/d)',ylab='SCFA (mM)',
            ylim=c(0,110),
             main=paste('RS frac ',RSf,' (cts inflow)'),
            cex.names=1.2,cex.axis=1.2,cex.lab=1.2,cex.main=1.4
            )
    legend('topleft',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.2,bty='n')

    barplot(t(scfa.array.all[,,'meals',RSf]),
            names.arg=CfromdietVec,
            col=scfa.cols,
            xlab='C inflow (g/d)',ylab='SCFA (mM)',
            ylim=c(0,110),
            main=paste('RS frac ',RSf,'d (meals)'),
            cex.names=1.2,cex.axis=1.2,cex.lab=1.2,cex.main=1.4
            )
    legend('topleft',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.2,bty='n')
}


dev.copy2eps(file='modelledDuncanExpt.eps')


#Duncan Expt SCFA
Duncan=rbind(M=c(65.09,17.67,20.28),
    HPMC=c(43.54,8.9,13.84),
    HPLC=c(35.50,4.36,10.81))
colnames(Duncan)=SCFAnames
DuncanNames=c('HPLC','HPMC','M')

#barplot(t(Duncan[c('HPLC','HPMC','M'),]),
#        col=scfa.cols,
#        main='Duncan et al.',
#        xlab='Diet',ylab='SCFA (mM)',
#        ylim=c(0,110),
#        cex.names=1.2,cex.axis=1.2,cex.lab=1.2,cex.main=1.4
#        )
#legend('topright',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.2,bty='n')

#dev.copy2eps(file='SCFAproteinExpt.eps')






#estimates for duncan data
DuncanDiet=rbind('M'=c(187,28),'HPMC'=c(95,12),'HPLC'=c(2.7,6.1))
colnames(DuncanDiet)=c('Starch','NSP')
DuncanC=NULL
for (i in 1:3){
    
    DuncanC=rbind(DuncanC,c(0,0.2,0.1)*DuncanDiet[i,'Starch']+0.75*DuncanDiet[i,'NSP'])
}
rownames(DuncanC)=rownames(DuncanDiet)
colnames(DuncanC)=c('lowC','highC','meanC')



dev.new(height=7,width=14)
par(mfrow=c(1,3),mar=c(6,6,2,1))
plot(c(0,60),c(0,30),type='n',xlab='C (g/d)',ylab='Butyrate (mM)',cex.lab=2,cex.axis=2)
ct=1
for (RSf in RSfvec){
    for (i in c('cts','meals')){

        mat=scfa.array.all[,,i,RSf]

        points(as.numeric(rownames(mat)),mat[,'Butyrate'],
               col=c('blue','red','blue','red')[ct],pch=c(15,17,15,17)[ct],
               cex=1.5)
        lines(as.numeric(rownames(mat)),mat[,'Butyrate'],
               col=c('blue','red','blue','red')[ct],lty=c(1,1,2,2)[ct],lwd=2,
               cex=1.5)

        ct=ct+1
    }
}
for (d in DuncanNames){
    lines(DuncanC[d,c('lowC','highC')],rep(Duncan[d,'Butyrate'],2),col='magenta',lwd=2)
}
lines(DuncanC[DuncanNames,'meanC'],Duncan[DuncanNames,'Butyrate'],lty=3,col='magenta',lwd=3)
points(DuncanC[DuncanNames,'meanC'],Duncan[DuncanNames,'Butyrate'],lty=3,col='magenta')

legend('topleft',col=c('blue','red','magenta'),pch=c(15,17,1),c('cts','meals','Duncan et al'),cex=2,bty='n',lwd=2)

legend('left',col='blue',lty=1:2,RSfvec,cex=1.7,bty='n',title='RS frac',lwd=2)

plot(c(0,60),c(0,40),type='n',xlab='C (g/d)',ylab='%Butyrate',cex.lab=2,cex.axis=2)
ct=1
for (RSf in RSfvec){
    for (i in c('cts','meals')){

        mat=scfa.array.all[,,i,RSf]
        tscfa1=rowSums(scfa.array.all[,,i,RSf])

        points(as.numeric(rownames(mat)),100*mat[,'Butyrate']/tscfa1,
               col=c('blue','red','blue','red')[ct],pch=c(15,17,15,17)[ct],
               cex=1.5)
        lines(as.numeric(rownames(mat)),100*mat[,'Butyrate']/tscfa1,
               col=c('blue','red','blue','red')[ct],lty=c(1,1,2,2)[ct],lwd=2,
               cex=1.5)

        ct=ct+1
    }
}
for (d in DuncanNames){
    lines(DuncanC[d,c('lowC','highC')],rep(100*Duncan[d,'Butyrate']/sum(Duncan[d,]),2),col='magenta',lwd=2)
}
lines(DuncanC[DuncanNames,'meanC'],100*Duncan[DuncanNames,'Butyrate']/rowSums(Duncan[DuncanNames,]),lty=3,col='magenta',lwd=3)
points(DuncanC[DuncanNames,'meanC'],100*Duncan[DuncanNames,'Butyrate']/rowSums(Duncan[DuncanNames,]),lty=3,col='magenta')


legend('topleft',col=c('blue','red','magenta'),pch=c(15,17,1),c('cts','meals','Duncan et al'),cex=2,bty='n',lwd=2)
legend('left',col='blue',lty=1:2,RSfvec,cex=1.7,bty='n',title='RS frac',lwd=2)

plot(c(0,60),c(0,120),type='n',xlab='C (g/d)',ylab='TSCFA (mM)',cex.lab=2,cex.axis=2)
ct=1
for (RSf in RSfvec){
    for (i in c('cts','meals')){

        mat=rowSums(scfa.array.all[,,i,RSf])

        points(CfromdietVec,mat,
               col=c('blue','red','blue','red')[ct],pch=c(15,17,15,17)[ct],
               cex=1.5)
        lines(CfromdietVec,mat,
               col=c('blue','red','blue','red')[ct],lty=c(1,1,2,2)[ct],lwd=2,
               cex=1.5)

        ct=ct+1
    }
}
#add in Duncan data
for (d in DuncanNames){
    lines(DuncanC[d,c('lowC','highC')],rep(sum(Duncan[d,]),2),col='magenta',lwd=2)
}
lines(DuncanC[DuncanNames,'meanC'],rowSums(Duncan[DuncanNames,]),lty=3,col='magenta',lwd=3)
points(DuncanC[DuncanNames,'meanC'],rowSums(Duncan[DuncanNames,]),lty=3,col='magenta')
      
legend('topleft',col=c('blue','red','magenta'),pch=c(15,17,1),c('cts','meals','Duncan et al'),cex=2,bty='n',lwd=2)
legend('bottomright',col='blue',lty=1:2,RSfvec,cex=1.7,bty='n',title='RS frac',lwd=2)


dev.copy2eps(file=paste0('ButyrateCProteinExpt.eps'))




    DuncBac=cbind(c(28.84,29.96,29.85),c(13.48,11.39,9.98),c(10.84,10.00,7.02),c(4.01,2.09,1.87),
        c(21.14,24.40,24.48),c(8.51,10.76,8.59),c(0.41,0.11,0.15),c(0.11,0.13,0.18),
        c(11.40,7.79,3.32),c(9.74,16.61,21.16),c(10.71,10.55,10.56)) 
    rownames(DuncBac)=c('M','HPMC','HPLC')
    colnames(DuncBac)=c('Bac303', 'Fprau645', 'Rfla729+Rbro730','Bif164', 'Erec482', 'Prop853', 'Lab158', 'Dsv698', 'Rrec584', 'Erec-Rrec', 'Total log count/g feces')


    DuncMat=cbind(
        B=DuncBac[,'Bac303'],
        NBSD=0.5*DuncBac[,'Rfla729+Rbro730']+ 0.3*(DuncBac[,'Erec482']-DuncBac[,'Rrec584']),
        NBFD=0.5*DuncBac[,'Rfla729+Rbro730']+ 0.3*(DuncBac[,'Erec482']-DuncBac[,'Rrec584']),
        LP=DuncBac[,'Bif164'],# + DuncBac[,'Ato291']
        BP1=DuncBac[,'Rrec584'],
        BP2=DuncBac[,'Fprau645'],
        PP=DuncBac[,'Prop853'],
        BP3=0.3*(DuncBac[,'Erec482'] - DuncBac[,'Rrec584']),
        A = 0.1*(DuncBac[,'Erec482']-DuncBac[,'Rrec584']),
        M =  0.001*DuncBac[,'Total log count/g feces']
    )

    maintenance.scfa=Duncan['M',]
    HPMC.scfa=Duncan['HPMC',]
    HPLC.scfa=Duncan['HPLC',]
    maintenance.bac=DuncMat['M',microbeNames.short]
    HPMC.bac=DuncMat['HPMC',microbeNames.short]
    HPLC.bac=DuncMat['HPLC',microbeNames.short]
