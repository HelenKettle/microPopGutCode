#after running 'makeCriteriaTable.R' load the validation table
rm(list=ls())
graphics.off()

library(ggplot2)

load('validationTable.Rdata')

boxnames=c('proximal','transverse','distal')

#reshape the data
TSCFA=Table[1:3,]
colnames(TSCFA)=rep(paste0(seq(1:3),'d'),2)
rownames(TSCFA)=c('prox','trans','dist')
PH=Table[4:6,]

mat=matrix(NA,ncol=5,nrow=18)
colnames(mat)=c('var','colon','TransitTime','inflow','value')

ct=1
for (m in c('meals','cts')){

    if (m =='cts'){
        tscfa=TSCFA[,1:3]
    }else{
        tscfa=TSCFA[,4:6]
    }
    rownames(tscfa)=boxnames
    colnames(tscfa)=paste0(seq(1:3),'d')

    for (box in boxnames){

        for (Tt in 1:3){
            
            mat[ct,]=c('TSCFA',box,Tt,m,tscfa[box,paste0(Tt,'d')])

            ct=ct+1
        }
    }
}


df=data.frame(mat)

print(df)

df$colon <- factor(df$colon, levels = boxnames)
df$value = as.numeric(mat[,'value'])

#dev.new()
ideal=data.frame(colon=boxnames,value=c(123,117,80))

d=ggplot(data=df,aes(x=TransitTime,y=value,color=inflow))+
    #geom_col(position="dodge") +
    geom_point(size=4) +
        facet_wrap(~colon) +
        geom_hline(data=ideal,aes(yintercept=value),color='deeppink',size=1.25) +
            labs(x='Transit time (d)',y='TSCFA (mM)')+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16,face="bold"),
                      strip.text.x = element_text(size = 20),
                      legend.title=element_text(size=20),
                      legend.text=element_text(size=14))
print(d)

dev.copy2eps(file='TSCFAvalidate.eps')

####-------------pH------------------------

PH=Table[4:6,]

mat.ph=matrix(NA,ncol=5,nrow=18)
colnames(mat.ph)=c('var','colon','TransitTime','inflow','value')

ct=1
for (m in c('meals','cts')){

    if (m =='cts'){
        ph=PH[,1:3]
    }else{
        ph=PH[,4:6]
    }
    rownames(ph)=boxnames
    colnames(ph)=paste0(seq(1:3),'d')

    for (box in boxnames){

        for (Tt in 1:3){
            
            mat.ph[ct,]=c('pH',box,Tt,m,ph[box,paste0(Tt,'d')])

            ct=ct+1
        }
    }
}


df.ph=data.frame(mat.ph)

print(df.ph)

df.ph$colon <- factor(df.ph$colon, levels = boxnames)
df.ph$value = as.numeric(mat.ph[,'value'])

dev.new()
ideal=data.frame(colon=boxnames,value=c(5.7,6.2,6.6))


d=ggplot(data=df.ph,aes(x=TransitTime,y=value,color=inflow))+
    geom_point(size=4) + facet_wrap(~colon) +
        geom_hline(data=ideal,aes(yintercept=value),size=1.25,color='deeppink') +
            labs(x='Transit time (d)',y='pH')+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16,face="bold"),
                          strip.text.x = element_text(size = 20),
                      legend.title=element_text(size=20),
                      legend.text=element_text(size=14))
print(d)

dev.copy2eps(file='PHvalidate.eps')

WF=Table[8,]
mat.w=matrix(NA,ncol=4,nrow=6)
colnames(mat.w)=c('var','TransitTime','meals','value')

ct=1
for (m in c('meals','cts')){

    if (m =='cts'){
        w=WF[1:3]
    }else{
        w=WF[4:6]
    }
    
    names(w)=paste0(seq(1:3),'d')
    
    for (Tt in 1:3){
        
        mat.w[ct,]=c('WF',Tt,m,w[paste0(Tt,'d')])
        
        ct=ct+1
    }
}

dev.new()
plot(as.numeric(mat.w[,'TransitTime']),100*as.numeric(mat.w[,'value']),type='n',xlab='Transit time (d)',ylab='% water remaining',cex.lab=1.5)

ii=mat.w[,'meals']=='meals'
points(as.numeric(mat.w[ii,'TransitTime']),100*as.numeric(mat.w[ii,'value']),col='blue',pch=15,cex=5,lwd=4)

ii=mat.w[,'meals']=='cts'
points(as.numeric(mat.w[ii,'TransitTime']),100*as.numeric(mat.w[ii,'value']),col='red',pch=21,cex=5,lwd=4)

legend('topright',legend=c('cts','meals'),col=c('red','blue'),pch=c(21,15),bty='n',cex=1.5)
