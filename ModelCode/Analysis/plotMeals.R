#plot the incoming meals for different scenarios

rm(list=ls())
graphics.off()

library(deSolve)

source('../../Rpackage/microPopGut/R/inflowFuncs.R')

numDays=2
time.step=1/(24*60)
RS.frac=0.78

C.from.diet=50
P.from.diet=10
W.from.diet=1100

RS.from.diet=RS.frac*C.from.diet
NSP.from.diet=(1-RS.frac)*C.from.diet

substrateMeans=c('P'=P.from.diet,'RS'=RS.from.diet,'NSP'=NSP.from.diet,'W'=W.from.diet)
substrateNames=names(substrateMeans)

fluc.subst.comp=TRUE
fluc.inflow=TRUE

Nseed=3

time=seq(0,numDays,time.step)
meals.array=array(NA,dim=c(length(time),length(substrateMeans),Nseed),dimnames=list(NULL,substrateNames,NULL))


scale=c(0.5*substrateMeans[c('P','RS','NSP')],0.1*substrateMeans['W'])

#plot Gamma distribution
x <- seq(-2, 150, length=10000)
cols=rainbow(length(scale))
names(cols)=substrateNames

subst.cols=c(P='red',RS='purple',NSP='black',W='blue')

dev.new()
par(mar=c(5,6,4,1))
ct=1
plot(range(x), c(0,0.08), type="n", xlab="inflow (g/d)",
  ylab="Density", main="a) gamma distributions",cex.lab=2,cex.axis=2,cex.main=2)
for (s in c('P','RS','NSP')){
    pdf=dgamma(x,shape=substrateMeans[s]/scale[s],scale=scale[s])
    lines(x,pdf, lwd=6, col=subst.cols[s],lty=ct)
    ct=ct+1
}
legend("topright", inset=.05, lwd=4, lty=1:(ct-1),col=subst.cols,legend=c('P','RS','NSP'),cex=1.5)
dev.copy2eps(file='inflowGamma.eps')


breakfast.start=7
lunch.start=13
dinner.start=19
meal.duration.h=0.5

time.step.h=24*time.step

for (seed in 1:Nseed){

    set.seed(seed)

    subst.mat=makeSubstrateMat(
        substrateNames,
        substrateMeans=substrateMeans,
        num.days=numDays,
        breakfast.start=breakfast.start,
        lunch.start=lunch.start,
        dinner.start=dinner.start,
        time.step.h=time.step*24,
        meal.duration.h=meal.duration.h,
        fluc.comp=fluc.subst.comp,
        time.to.reach.colon.h=7,
        gamma.mag=1,
        RS.frac=RS.frac,
        showPlot=FALSE,
        waterName='H2O',
        waterMean=substrateMeans['W']
    )

    
    meals.array[,,seed]=subst.mat[,2:5]

    time=subst.mat[,1]

}

seed.cols=c('orange','cyan','darkgreen')

dev.new()
par(mfrow=c(2,2),mar=c(4,6,4,1))

for (sn in substrateNames){
    if (sn=='P'){
        main.txt=paste0('b) Inflow for ',Nseed,' random examples')
    }else{
        main.txt=''
    }
    matplot(time,meals.array[,sn,],
            type='l',ylab=paste0(sn,' inflow (g/d)'),xlab='time (d)',
            main=main.txt,lwd=4,lty=1:3,
            cex.lab=1.5,cex.axis=1.2,col=seed.cols,cex.main=1.2)
    legend('topleft',legend=paste(1:Nseed),lwd=4,lty=1:3,col=seed.cols,bty='n',cex=1,title='seed')
}
dev.copy2eps(file='inflowTS.eps')


#plot composition
dev.new(height=7,width=9)
par(mfrow=c(2,Nseed),mar=c(5,5,5,2))
Nmeals=numDays*3
for (i in 1:2){
    if (i==1){snames=c('P','NSP','RS');cols=subst.cols[c('P','NSP','RS')];ymax=150}
    if (i==2){snames='W';subst.cols['W'];ymax=1600}
    Ns=length(snames)
    for (seed in 1:Nseed){
        set.seed(seed)
        meals.only=matrix(NA,nrow=Nmeals,ncol=Ns)
        colnames(meals.only)=snames
        for (meal in 1:Nmeals){
            for (s in snames){
                if (substrateMeans[s]>0){
                    meals.only[meal,s]=rgamma(1,shape=substrateMeans[s]/scale[s],scale=scale[s])
                }else{
                    meals.only[meal,s]=0
                }
            }
        }
    #rescale to maintain the average
        meals.scaled=meals.only*NA
        for (s in snames){
            if (substrateMeans[s]>0){
                meals.scaled[,s]=substrateMeans[s]*meals.only[,s]/mean(meals.only[,s])
            }else{
                meals.scaled[,s]=0
            }
            #print(paste(s,':SD=',sqrt(var(meals.scaled[,s]))))
        }

        if (seed==1 & i==1){
            main.txt=paste0('c) Example composition of ',Nmeals,' meals')
        }else{
            main.txt=' '
        }
        barplot(t(meals.scaled[,snames]),col=subst.cols[snames],main=main.txt,
                ylim=c(0,ymax),legend.text=snames,names.arg=1:Nmeals,
                xlab='meals',ylab='inflow (g/d)',
                cex.lab=1.5,cex.main=1.4,cex.axis=1.5,cex.names=1.5)
        
        mtext(paste('seed =',seed))

        #print(colMeans(meals.only))
        #print(colMeans(meals.scaled))
    }
}

dev.copy2eps(file='inflowComp.eps')


