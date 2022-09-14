
#creates blank plot
blankPlot=function(withLegend=FALSE,...){
    plot(c(0,1),c(0,1),type='n',xlab='',ylab='',xaxt='n',bty='n',yaxt='n')
    if (withLegend){
        legend('top',col=cols,legend=text,lty=ltype,lwd=lwd)
    }
}
