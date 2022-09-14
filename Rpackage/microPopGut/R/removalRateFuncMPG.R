#' Removal Rate Function for microPopGut
#' 
#' Return the rate of removal of any state variable from the system
#' This is called in the ODE derivs func
#' 
#' @aliases removalRateFunc
#'
#' @param varName (string) Name of state variable of interest (this is group name
#' or a resource name - NOT a strain name)
#' @param varValue (scalar) value of state variable of interest
#' @param stateVarValues (named vector) values of all state variables
#' @param time (scalar) time
#' @param washOut (named vector) of wash out rates (per unit time) of groups
#' and resources (specified in SysInfo files)
#' @param parms List containing all system parameters
#' @return (scalar) rate of removal (quantity per unit time) for the state variable varName
#' @export
#' 
removalRateFuncMPG=function (varName, varValue, stateVarValues, time, washOut, parms){

    if (varValue <= 0) {
        v = 0
    }else{
        v = washOut[varName] * varValue
    }

    if (parms$myPars$BMpd>0){
        if (parms$myPars[['compartment']]==parms$myPars$nbox){
            v=0
            GMT.h=(24*time)%%24
            start.times=parms$myPars$start.BM.time[[parms$myPars$BMpd]]
            for (b in 1:parms$myPars$BMpd){
                st=start.times[b]
                if (GMT.h>st & GMT.h<=(st + parms$myPars$BM.duration.h)){
                    v=parms$myPars$V.bm*varValue
                }
        
            }#b
        }

    }

    
    v.abs=0
    #SCFA absorption
    if (varName%in%parms$myPars$SCFAnames){
        v.abs=parms$myPars$abs.scfa*varValue
    }
    #water absorption
    if (varName==parms$myPars$waterName){
        v.abs=parms$myPars$abs.water*varValue
    }

    return(v+v.abs)
}

