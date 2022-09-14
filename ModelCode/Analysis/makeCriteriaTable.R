rm(list=ls())
graphics.off

library(xtable)
#library(brewer)

source('plotFile.R')

#folderName='/home/helen/MicroPop/Applications/Human/Aug21/Results/Experiments/New/Protein/'
folderName='/home/helen/MicroPop/Applications/Human/Aug21/Results/Experiments/New/transitTime/'


Ttvec=c('1','2','3')
Fvec=c('TRUE','FALSE')
Table=matrix(NA,ncol=length(Ttvec)*2,nrow=8)
STable=rep(NA,length(Ttvec)*2)
ct=1
for (f in 1:length(Fvec)){

    fluc.inflow=Fvec[f]
    
    if (fluc.inflow){Nseed=4}else{Nseed=1}
    
    for (it in 1:length(Ttvec)){
    
        Tt=Ttvec[it]

        col=matrix(NA,nrow=8,ncol=Nseed)
        A=rep(NA,Nseed)
        B=rep(NA,Nseed)
        P=rep(NA,Nseed)

        for (seed in 1:Nseed){
            
            file=paste0('numDays28seed',seed,'fix.pHFALSEfluc.inflow',fluc.inflow,'fluc.subst.comp',fluc.inflow,'numStrains1BMpd0C.from.diet50P.from.diet10W.from.diet1100transitTime',Tt,'runNumber1RS.frac0.78LPc7gamma.mag1.Rdata')
            
            criteria=plotFile(file,folderName,plot=FALSE)
        
            col[,seed]=c(criteria$TSCFA,criteria$pH,
                   fecalOutput=criteria$fecalOutput,W.frac=criteria$w.frac)
            A[seed]=criteria$SCFA['dist','Acetate']
            B[seed]=criteria$SCFA['dist','Butyrate']
            P[seed]=criteria$SCFA['dist','Propionate']
        }
        
        av1=rowMeans(col,na.rm=TRUE)

        avA=mean(A,na.rm=TRUE)
        avB=mean(B,na.rm=TRUE)
        avP=mean(P,na.rm=TRUE)
        
        Table[,ct]=av1
        
        STable[ct]=paste0(round(avA),':',round(avB),':',round(avP))
        ct=ct+1
        
    }
}


print(Table)

xtable(Table)
TSCFAmat=Table[1:3,]
colnames(TSCFAmat)=c('Meals.1d','Meals.2d','Meals.3d','Cts.1d','Cts.2d','Cts.3d')
rownames(TSCFAmat)=c('TSCFA.prox','TSCFA.trans','TSCFA.dist')

dev.new()
par(mfrow=c(1,3),mar=c(4,4,2,2))
plot(c(1,3),c(60,130),type='n',xlab='Transit time (d)',ylab='TSCFA (mM)',main='Prox')
points(TSCFAmat['TSCFA.prox',1:3],col='red')
points(TSCFAmat['TSCFA.prox',4:6],col='blue')

plot(c(1,3),c(60,130),type='n',xlab='Transit time (d)',ylab='TSCFA (mM)',main='Trans')
points(TSCFAmat['TSCFA.trans',1:3],col='red')
points(TSCFAmat['TSCFA.trans',4:6],col='blue')

plot(c(1,3),c(60,130),type='n',xlab='Transit time (d)',ylab='TSCFA (mM)',main='Dist')
points(TSCFAmat['TSCFA.dist',1:3],col='red')
points(TSCFAmat['TSCFA.dist',4:6],col='blue')

