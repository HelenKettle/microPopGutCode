#plot Curves for fitting pH, KpH, TSCFA etc
#graphics.off()

pH.colon=c(5.7,6.2,6.6)
TSCFA.colon=c(123,117,80)
TSCFAmin=80

dev.new(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,5,2,1))
S=seq(0,200)
plot(range(S),c(4,9),type='n',xlab='TSCFA (mM)',ylab='pH',
     cex.axis=1.5,cex.lab=1.8,main='a) relating pH to TSCFA')
#pH=6.6*(1-)^S
#lines(S,pH,col='red')
#lines(S,11.4*(1-(S/123)^2),col=2)
#lines(S,5.7+(1-((S-80)/43)^2),col=2)
points(TSCFA.colon,pH.colon,lwd=3)
#lines(TSCFA.colon,pH.colon)

lm1=lm(pH.colon~TSCFA.colon)

ph=lm1$coef[1]+lm1$coef[2]*S

#lines(S,ph)

pH=rep(NA,length(S))
for (i in 1:length(S)){
    TSCFA=S[i]
    if (TSCFA<80){
        pH[i]=6.65
    }else if (TSCFA>140){
        pH[i]=4
    }else{
        pH[i] = 6.65-((TSCFA-80)/43)^5.3
    }
}

pH=8.02-0.0174*S
pH[pH<5]=5
pH[pH>8]=8

#lines(S,pH,col='green')
#lines(S,8-0.5*S/40,col='magenta')
lines(S,pH,col='blue',lwd=3)


#plot pH limit function
plot(c(5,8),c(0,1),type='n',xlab='pH',ylab='scaling on microbial growth',
     cex.lab=1.8,cex.axis=1.5,main='b) microbial pH tolerance')
points(c(0,5.5,6.5,7,7.5,8.5),c(0,0,1,1,0,0))
lines(c(0,5.5,6.5,7,7.5,8.5),c(0,0,1,1,0,0),lwd=4,col='blue')

dev.copy2eps(file='pHTSCFAandTol.eps')

dev.copy2pdf(file='pHTSCFAandTol.pdf')

more=0
if (more){
avec=NULL
bvec=NULL
evec=NULL
for (a in seq(5,8,0.01)){
    for(b in seq(1,6,0.1)){
        p=a-((S-80)/43)^b
        error=sum((pH.colon-approx(S,p,TSCFA.colon)$y)^2)
        avec=c(avec,a)
        bvec=c(bvec,b)
        evec=c(evec,error)
    }
}

besta=avec[evec==min(evec)]
bestb=bvec[evec==min(evec)]

print(paste('a=',besta))
print(paste('b=',bestb))

lines(S,besta-((S-80)/43)^bestb,col=2)
dev.copy2eps(file='pHTSCFAfunction.eps')

dev.new()
TSCFA=seq(80,140)
pH=besta-((TSCFA-80)/43)^bestb
plot(TSCFA,pH)

TSCFA=seq(0,140)
pH=8-0.5*TSCFA/40

lines(TSCFA,pH,lwd=3)
points(TSCFA.colon,pH.colon,col='red')
}
