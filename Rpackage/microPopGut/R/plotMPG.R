#' plotMPG
#'
#' Basic plot of the model results
#'
#' @param m.out is the output of microPopGut(). It is a list containing the following objects: solution, params, myPars, etc
#' @export


plotMPG=function(m.out){

    
    molarMass=m.out$params[[1]]$molarMass
    SCFAnames=m.out$myPars[[1]]$SCFAnames
    waterName=m.out$myPars[[1]]$waterName
    allStrainNames=m.out$params[[1]]$allStrainNames
    microbeNames=m.out$params[[1]]$microbeNames

    m.col=rainbow(length(microbeNames))
    names(m.col)=microbeNames
    
    gNames.short=m.out$myPars[[1]]$microbeNames.short[m.out$params[[1]]$microbeNames]

    numStrains=m.out$params[[1]]$numStrains

    #convert numStrains to a vector if it isn't one

    if (length(numStrains)==1){
        numStrains=rep(numStrains,length(microbeNames))
        names(numStrains)=microbeNames
    }
        
    if (any(numStrains>1)){
        strainNames.short = NULL
        groupName = NULL
        for (g in microbeNames){
            strainNames.short=c(strainNames.short,paste0(gNames.short[g],'.',seq(1,numStrains[g])))
            groupName=c(groupName,rep(g,numStrains[g]))
        }
    }else{
        strainNames.short=gNames.short
        groupName=microbeNames
    }
    
    
    time=m.out$solution[[1]][,1]

    #dev.new()
    par(mfcol=c(3,3),mar=c(4,4,2,0.5))

    SCFA=list()
    microbes=list()
    pH=list()

    
    for (box in 1:3){
        
        SCFA[[box]]=computeSCFA(m.out$solution[[box]],SCFAnames,molarMass,waterName)
        
        microbes[[box]]=m.out$solution[[box]][,allStrainNames]

        pH[[box]]=m.out$solution[[box]][,'pH']

    }

    scfa.names=colnames(SCFA[[1]])
    
    for (box in 1:3){

        matplot(time,SCFA[[box]],type='l',col=1:3,lwd=2,
                main=paste(m.out$myPars[[box]]$compartmentNames[[box]]),
                ylab='SCFA (mM)',xlab='time (d)',
                ylim=c(0,1.2*max(unlist(SCFA))),
                cex.lab=1.2,cex.axis=1.2)
        legend('topleft',scfa.names,col=1:3,lwd=2,lty=1:3,bty='n')
        
        matplot(time,microbes[[box]],type='n',col=m.col,lwd=2,
                main=paste(m.out$myPars[[box]]$compartmentNames[[box]]),
                ylab='Microbes (g)',xlab='time (d)',
                ylim=c(0,1.2*max(unlist(microbes))),
                cex.lab=1.2,cex.axis=1.2)

        if (length(microbeNames)==1){
            lines(time,microbes[[box]],col=m.col[microbeNames])
        }else{
            for (g in microbeNames){
                matlines(time,microbes[[box]][,allStrainNames[groupName==g]],col=m.col[g])
            }
        }

        str.col=NULL
        str.lty=NULL
        for (g in microbeNames){
            str.col=c(str.col,rep(m.col[g],numStrains[g]))
            str.lty=c(str.lty,seq(1,numStrains[g]))
        }
        
        legend('topleft',strainNames.short,col=str.col,lwd=2,lty=str.lty,bty='n')

        plot(time,pH[[box]],type='l',xlab='time (d)',ylab='pH',lwd=2,
             ylim=c(min(unlist(pH)),max(unlist(pH))),
             main=paste(m.out$myPars[[box]]$compartmentNames[[box]]),
             cex.lab=1.2,cex.axis=1.2)

        
    }
    
}
