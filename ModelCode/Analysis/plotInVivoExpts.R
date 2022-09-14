rm(list=ls())
graphics.off()

library('RColorBrewer')

scfa.cols=brewer.pal(3,'Set1')

SCFAnames=c('Acetate','Butyrate','Propionate')

Duncan=rbind(M=c(65.09,17.67,20.28),
    HPMC=c(43.54,8.9,13.84),
    HPLC=c(35.50,4.36,10.81))
colnames(Duncan)=SCFAnames

DuncanNames=c('M','HPMC','HPLC')

#Walker SCFA data is from an email from Harry of some analysis by Grietje
scfa.walker=rbind(
    M=c(0.24,77.61,26.06,1.68,23.92,1.31,2.42,0.45,3.47,4.21,0.08,141.47),
    NSP=c(0.62,81.74,22.37,1.26,24.00,0.98,1.98,0.32,2.24,8.64,0.11,144.26),
    RS=c(0.75,54.61,15.45,1.65,14.82,1.33,1.94,0.86,3.70,1.27,0.16,96.53),
WL=c(0.11,59.78,17.07,1.78,12.97,1.42,2.01,0.54,3.85,0.90,0.08,100.49)
)
colnames(scfa.walker)=c('Formate','Acetate','Propionate','Isobutyrate','Butyrate','Isovalerate','Valerate','Cap','Lactate','Succinate','Fum','SCFAtot')

expts=c('M','NSP','RS','WL')
NSPI=c('M'=28,'NSP'=42,'RS'=13,'WL'=22)
RSI=c('M'=230,'NSP'=138,'RS'=275,'WL'=110)
CI=c('M'=427,'NSP'=427,'RS'=434,'WL'=201)


#RS.frac=RSI/(RSI+NSPI)
RS.frac=RSI/CI
exptsNew=expts[order(RS.frac)]

mat1=scfa.walker[exptsNew,SCFAnames]


dev.new(width=10,height=7)
par(mfrow=c(1,2),mar=c(5,5,2,1))
barplot(t(mat1),col=scfa.cols,names.arg=exptsNew,
        main='Walker et al. 2011',
        xlab='Diet',
        ylim=c(0,130),
        ylab='SCFA (mM)',
        cex.lab=1.5,cex.names=1.5,cex.main=1.5,cex.axis=1.5
        )
legend('topright',pch=15, col=scfa.cols,SCFAnames,cex=1.6,bty='n')


DuncanNamesNew=c('HPLC','HPMC','M')

barplot(t(Duncan[DuncanNamesNew,]),col=scfa.cols,names.arg=DuncanNamesNew,
        main='Duncan et al. 2007',
        xlab='Diet',
        ylim=c(0,130),
        ylab='SCFA (mM)',
        cex.lab=1.5,cex.names=1.5,cex.main=1.5,cex.axis=1.5
        )
legend('top',pch=15, col=scfa.cols,SCFAnames,cex=1.6,bty='n')

dev.copy2eps(file='InVivoExptsNew.eps')
