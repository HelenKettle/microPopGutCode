verification=function(m.out,start.av,fin.av){

    #' verification
    #'
    #' Prints a summary of the model results averaged over the time period start.av to fin.av
    #'
    #' @param m.out is the output of microPopGut(). It is a list containing the following objects: solution, params, myPars, etc
    #' @param start.av if the time to start averaging over
    #' @param fin.av if the time to finish averaging over
    #' @export
  

    molarMass=m.out$params[[1]]$molarMass
    SCFAnames=m.out$myPars[[1]]$SCFAnames
    waterName=m.out$myPars[[1]]$waterName
    compartNames=m.out$myPars[[1]]$compartmentNames
    time=m.out$solution[[1]][,'time']

    microbeNames=m.out$params[[1]]$allStrainNames
    
    microbes=matrix(NA,nrow=length(microbeNames),ncol=3,
        dimnames=list(microbeNames,compartNames))

    av.indices=seq(1,length(time))[time>=start.av & time<fin.av]

    SCFA=NULL
    for (box in 1:3){

        scfa=computeSCFA(m.out$solution[[box]],SCFAnames,molarMass,waterName)

        if (length(av.indices)>1){
            av.scfa=c(sum(colMeans(scfa[av.indices,])),colMeans(scfa[av.indices,]))
        }else{
            av.scfa=c(sum(scfa[av.indices,]),scfa[av.indices,])
        }

        SCFA=cbind(SCFA,av.scfa)
        
    }

    rownames(SCFA)=paste(c('TSCFA',names(av.scfa[2:length(av.scfa)])),'(mM)')
    colnames(SCFA)=compartNames
    
    print(round(SCFA,1))

    Vw=m.out$params[[3]]$Smats$washOut[waterName]
    Wfrac=mean(m.out$solution[[3]][av.indices,waterName])*Vw/m.out$myPars[[1]]$W.from.diet
    print(paste('Fraction of the incoming water that leaves colon is',round(Wfrac,2)))

    Vm=m.out$params[[3]]$Smats$washOut[m.out$params[[1]]$microbeNames]
    if (length(av.indices)>1){
        if (length(microbeNames)>1){
            fecal.output=mean(rowSums(m.out$solution[[3]][av.indices,microbeNames]))*mean(Vm)
        }else{
            fecal.output=mean(m.out$solution[[3]][av.indices,microbeNames])*mean(Vm)
        }
    }else{
        if (length(microbeNames)>1){
            fecal.output=sum(m.out$solution[[3]][av.indices,microbeNames])*mean(Vm)
        }else{
            fecal.output=m.out$solution[[3]][av.indices,microbeNames]*mean(Vm)
        }
    }
    print(paste('fecal microbe output rate is',round(fecal.output,2),'g/d'))

}
