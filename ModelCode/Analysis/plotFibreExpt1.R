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

lead.var='RS.frac'

folderName=paste0(resFolder,'Fibre/')

all.files=list.files(folderName)
Lf=length(all.files)

print(paste0('there are ',Lf,' files'))

ft=makeFileTable(all.files,batch.par.names)

#CfromdietVec=
#PfromdietVec=c(M=10,RS=10,NSP=10)

#numDaysVec=c(M=14,RS=21,NSP=21)
#numDaysVec=c(M=28,RS=28,NSP=28)
#RSfracVec=c(M=0.2,RS=0.68,NSP=0.074)
#RSfracVec=c(M=0.78,RS=0.95,NSP=0.05)


#RSfracVec['NSP']=0.045
#CfromdietVec['NSP']=44
Ttvec=paste(c(1,3))
#diets=c('M','RS','NSP')

RS.fracVec=c(0.1,0.5,0.78,0.9)
diets=paste(RS.fracVec)

bac.array.all=array(NA,dim=c(length(diets),10,2,length(Ttvec)),dimnames=list(diets,microbeNames.short,c('cts','meals'),paste(Ttvec)))
scfa.array.all=array(NA,dim=c(length(diets),3,2,length(Ttvec)),dimnames=list(diets,SCFAnames,c('cts','meals'),paste(Ttvec)))

for (Tt in Ttvec){


    all.bac=array(NA,dim=c(length(diets),10,2),
        dimnames=list(diets,microbeNames.short,c('cts','av.seed')))

    all.scfa=array(NA,dim=c(length(diets),3,2),
        dimnames=list(diets,SCFAnames,c('cts','av.seed')))

    for (diet in diets){
    
        numDays=28#numDaysVec[diet]

#find the files I want
        fix.var=c('numDays'=paste(numDays),
            transitTime=Tt,
            RS.frac=diet
                                        #C.from.diet=unname(CfromdietVec[diet]),
                                        #P.from.diet=unname(PfromdietVec[diet])
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
            dev.new(height=length(diets)*3,width=(num.other.var.vals+1)*3)
            par(mfrow=c(length(diets),num.other.var.vals+1),c(4,5,2,0.5))
            dev.new(height=length(diets)*3,width=(num.other.var.vals+1)*3)
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
                main.txt=paste0('RS frac = ',inf[lead.var],', meals (seed=',inf['seed'],')')
            }else{
                main.txt=paste0('RS frac = ',inf[lead.var],', continuous inflow')
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

    bac.array.all[,,,paste(Tt)]=all.bac
    scfa.array.all[,,,paste(Tt)]=all.scfa

}#Tt

#plot mean SCFA
dev.new(height=7,width=10)
par(mfrow=c(2,length(Ttvec)*2),mar=c(5,5,4,1))

for (Tt in Ttvec){
    barplot(t(bac.array.all[,,'cts',Tt]),
            col=bac.cols,
            xlab='RS fraction',ylab='Biomass (g)',
            ylim=c(0,35),
             main=paste('transit time',Tt,'d (cts inflow)'),
            cex.names=1.4,cex.axis=1.4,cex.lab=1.4
            )
    if (Tt==Ttvec[1]){
        legend('topleft',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.1,bty='n')
    }
    
    barplot(t(bac.array.all[,,'meals',Tt]),
            col=bac.cols,
            xlab='RS fraction',ylab='Biomass (g)',
            ylim=c(0,35),
            main=paste('transit time',Tt,'d (meals)'),
            cex.names=1.4,cex.axis=1.4,cex.lab=1.4
            )
    if (Tt==Ttvec[1]){
        legend('topleft',col=bac.cols,legend=microbeNames.short,pch=15,cex=1.1,bty='n')
    }

}

for (Tt in Ttvec){
    barplot(t(scfa.array.all[,,'cts',Tt]),
            col=scfa.cols,
            xlab='RS fraction',ylab='SCFA (mM)',
            ylim=c(0,120),
             main=paste('transit time',Tt,'d (cts inflow)'),
            cex.names=1.4,cex.axis=1.4,cex.lab=1.4
            )
    legend('topleft',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.1,bty='n')
    
    
    barplot(t(scfa.array.all[,,'meals',Tt]),
            col=scfa.cols,
            xlab='RS fraction',ylab='SCFA (mM)',
            ylim=c(0,120),
            main=paste('transit time ',Tt,'d (meals)'),
            cex.names=1.4,cex.axis=1.4,cex.lab=1.4
            )
    legend('topleft',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.1,bty='n')
}

dev.copy2eps(file='modelledRSfrac.eps')



dev.new(width=10,height=7)
par(mfrow=c(1,length(Ttvec)+1),mar=c(5,5,4,1))

#pick high RS and Tt=3
#low RS and Tt=1
#to represent how fibre changes Tt
for (i in c('cts','meals')){
    matM=rbind(scfa.array.all[c('0.1','0.5'),,i,'1'],
        scfa.array.all[c('0.78','0.9'),,i,'3'])
    barplot(t(matM),
            col=scfa.cols,
            xlab='RS fraction',ylab='SCFA (mM)',
            ylim=c(0,130),
            main=paste(i,'inflow with variable Tt'),
            cex.lab=1.3,cex.names=1.3,cex.main=1.5,cex.axis=1.3
            )
    if (i=='cts'){
        legend('topright',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.5)
    }
}
scfa.walker=rbind(
    M=c(0.24,77.61,26.06,1.68,23.92,1.31,2.42,0.45,3.47,4.21,0.08,141.47),
    NSP=c(0.62,81.74,22.37,1.26,24.00,0.98,1.98,0.32,2.24,8.64,0.11,144.26),
    RS=c(0.75,54.61,15.45,1.65,14.82,1.33,1.94,0.86,3.70,1.27,0.16,96.53))
colnames(scfa.walker)=c('Formate','Acetate','Propionate','Isobutyrate','Butyrate','Isovalerate','Valerate','Cap','Lactate','Succinate','Fum','SCFAtot')

expts=c('NSP','M','RS')
mat1=scfa.walker[expts,SCFAnames]
barplot(t(mat1),col=scfa.cols,names.arg=expts,
        main='Walker',
        xlab='Diet',
        ylim=c(0,130),
        ylab='SCFA (mM)',
        cex.lab=1.3,cex.names=1.3,cex.main=1.5,cex.axis=1.3
        )
dev.copy2eps(file='SCFA_RSfracData.eps')



#plot bacterial shifts
mdiet='0.78'
for (diet in c('0.1','0.9')){
    dev.new(height=9,width=7)
    par(mfrow=c(length(Ttvec),2),mar=c(5,5,2,1))
    for (Tt in paste(Ttvec)){
        for (i in c('cts','meals')){
            mat1=log(bac.array.all[diet,,i,Tt]/bac.array.all[2,,i,Tt])
            barplot(mat1,horiz=TRUE,col=bac.cols,names.arg=microbeNames.short,
                    las=1,
                    main=paste0('RS frac ',mdiet,' -> ',diet,' (Tt=',Tt,', ',i,')'),
                    xlab=paste0('log(M',diet,'/M',mdiet,')'),
                    xlim=c(-15,10)
                    )
        }
    }
    dev.copy2eps(file=paste0('MicrobeSwingRSfrac',diet,'.eps'))
}




more=0
if (more){        
#PCR.names=c('Bac','Ros','Fprau','Rum','Osc','Bif','Met')
#convert to microPop groups
#note: not sure if 'Rum' is NBFD or NBSD
#PCR.names=c('Bacteroides','ButyrateProducers1','ButyrateProducers2','NoButyFibreDeg','NoButyStarchDeg','LactateProducers','Methogens')
PCR.names=c('B','BP1','BP2','NBFD','NBSD','LP','M')
PCR.M=c(27.8,7.3,11.2,6.5,0.74,1.9,0.080)
PCR.NSP=c(25.7,6.5,14.4,3.8,0.77,1.8,0.034)
PCR.RS=c(20.2,10.1,12.1,17.0,2.0,2.4,0.124)

PCR.walker=rbind(PCR.M,PCR.RS,PCR.NSP)
colnames(PCR.walker)=PCR.names
rownames(PCR.walker)=c('M','RS','NSP')

#if Rum is NBSD
PCR.names=c('B','BP1','BP2','NBSD','LP','M')
PCR.M=c(27.8,7.3,11.2,6.5+0.74,1.9,0.080)
PCR.NSP=c(25.7,6.5,14.4,3.8+0.77,1.8,0.034)
PCR.RS=c(20.2,10.1,12.1,17.0+2.0,2.4,0.124)

PCR.walker=rbind(PCR.M,PCR.RS,PCR.NSP)
colnames(PCR.walker)=PCR.names
rownames(PCR.walker)=c('M','RS','NSP')

walker.mat=matrix(0,ncol=10,nrow=3)
colnames(walker.mat)=microbeNames.short
rownames(walker.mat)=c('M','RS','NSP')

walker.mat[,PCR.names]=PCR.walker


#SCFA measurements (means)

scfa.walker=rbind(
    M=c(0.24,77.61,26.06,1.68,23.92,1.31,2.42,0.45,3.47,4.21,0.08,141.47),
    NSP=c(0.62,81.74,22.37,1.26,24.00,0.98,1.98,0.32,2.24,8.64,0.11,144.26),
    RS=c(0.75,54.61,15.45,1.65,14.82,1.33,1.94,0.86,3.70,1.27,0.16,96.53))
colnames(scfa.walker)=c('Formate','Acetate','Propionate','Isobutyrate','Butyrate','Isovalerate','Valerate','Cap','Lactate','Succinate','Fum','SCFAtot')

dev.new(height=7,width=9)
par(mfrow=c(1,3),mar=c(5,5,2,1))
for (i in c('cts','meals')){
    mat1=scfa.walker[c('M','RS','NSP'),SCFAnames]
    barplot(t(mat1),col=scfa.cols,names.arg=c('M','RS','NSP'),
            main='Walker',
            xlab='Diet',
            ylim=c(0,150),
            ylab='SCFA (mM)',
            cex.lab=1.3,cex.names=1.2
                    )
legend('topright',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.1)
}

mat1=scfa.walker[c('M','RS','NSP'),c('Acetate','Butyrate','Propionate')]
barplot(t(mat1),col=scfa.cols,names.arg=c('M','RS','NSP'),
        main=paste0('Walker'),
        ylim=c(0,150),
        ylab='SCFA (mM)',
        cex.lab=1.3,cex.names=1.2
        )

legend('topright',col=scfa.cols,legend=SCFAnames,pch=15,cex=1.1)
}
