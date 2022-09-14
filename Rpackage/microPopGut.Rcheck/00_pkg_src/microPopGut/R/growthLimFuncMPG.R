#' growth rate limitation function for microPopGut
#' 
#' Returns the value of growthLim (must lie in interval \[0,1] i.e. unitless) of strainName on varName which is used to scale the maximum growth rate
#' Contains two options - one for essential resources and one for substitutable resources (based on Ballyk and Wolkowicz, 1993)
#'
#' @aliases growthLimFunc

#' @param strainName Name of the strain that is being looped through in the ODE solver
#' @param groupName Name of microbial group that is being looped through in the ODE solver
#' @param pathName Name of metabolic path (e.g. path1) that is being looped through in the ODE solver
#' @param varName (string) Name of variable (resource) of interest 
#' @param resourceValues State vector of resources  (with names)
#' @param allSubType Vector of strings (with names corresponding to the resourceNames) which describes the type of each resource ('Rtype') - Rtypes are S (substitutable resource), Se (essential resource), Sb (booster resource), Sm (microbial resource), P (product) and Pb (biomass product)
#' @param strainHalfSat Vector (with names corresponding to the resourceNames) of half-saturation constants for the given strain. If resource is not a substrate for the given strain, the value is NA 
#' @param stateVarValues State vector (resources and microbes) (with names)
#' @param parms list of parameter values
#' 
#' @return scalar giving limitation on growth rate - must be >=0 and <=1
#' @export
#' 
growthLimFuncMPG = function (strainName, groupName, pathName,
    varName, resourceValues, allSubType, strainHalfSat, stateVarValues, parms){


    v=1

    #added this, in case water is not a resource of the microbe
    if (!parms$myPars$waterName%in%names(allSubType)){
        allSubType=c(allSubType,"X")
        names(allSubType)=c(names(allSubType[1:(length(allSubType)-1)]),parms$myPars$waterName)
    }
    
    
    if (resourceValues[varName] <= 0) {

        v = 0
        
    }else {
        
        Wvol=stateVarValues[parms$myPars$waterName]/1000

        
        if (allSubType[varName] == "Sb" | allSubType[varName] == "Se") {
            v = resourceValues[varName]/(resourceValues[varName] + 
                                             strainHalfSat[varName]*Wvol)
        }
        
        if (allSubType[varName] == "S" | allSubType[varName] == "Sm") {

            S.res=names(allSubType[allSubType =='S'])
            Sm.res=names(allSubType[allSubType =='Sm'])

            Sr=NULL
            if (length(S.res)>=1){
                Sr=c(Sr,S.res)
            }
            if (length(Sm.res)>=1){
                Sr=c(Sr,Sm.res)
            }

            
            v = resourceValues[varName]/(Wvol*strainHalfSat[varName] * 
                (1 + sum(resourceValues[Sr]/
                 (Wvol*strainHalfSat[Sr]))))

            
        }
    }

    return(max(v, 0))
}

