#' entry Rate Function for microPopGut
#' 
#' Return the rate of entry to the system for any state variable
#'
#' @aliases entryRateFunc
#'
#' @param varName (string) Name of state variable of interest (resource name or strain name)
#' @param varValue (scalar) value of state variable of interest
#' @param stateVarValues (named vector) values of all state variables
#' @param time (scalar) time
#' @param inflowRate (named vector) on inflow rates (specified in SysInfo files)
#' @param parms List containing all system parameters
#' @return (scalar) rate of entry (quantity per unit time) for any state variable
#'
#' @export
entryRateFuncMPG=function(varName,varValue,stateVarValues,time,inflowRate,parms){


  #output is the resource (or microbial strain mass) per unit time

    if (as.numeric(parms$myPars[['compartment']])==1){

        gname = microPop::getGroupName(varName, parms$microbeNames)

        if (gname %in% parms$microbeNames) {
            if (length(parms$numStrains)==1){
                v = inflowRate[gname]/parms$numStrains
            }else{
                v = inflowRate[gname]/parms$numStrains[gname]
            }
            
        } else {
        
            if (parms$myPars$fluc.inflow & varName%in%c(parms$myPars$waterName,parms$myPars$substrateNames)){
                #use time series for subsrate inflow
                v=stats::approx(parms$myPars$subst.inflow[,1],parms$myPars$subst.inflow[,varName],time,rule=2)$y
                
            }else{
                v=inflowRate[varName]
            }
        }
        
    }else if (parms$myPars[['compartment']]>1){
   #now the entry rate is the output from the previous compartment at the given time step
        index=stats::approx(parms$myPars$inflow.mat[,'time'],seq(1,nrow(parms$myPars$inflow.mat)),time,method='constant',rule=2)$y
        v=parms$myPars$inflow.mat[index,varName]
        
    }

    if (varName%in%c(parms$myPars$carbs,'Protein')){

        Wvol=stateVarValues[parms$myPars$waterName]/1000


        if (!is.null(parms$myPars$mucin.degrader)){


            if (length(parms$numStrains)==1){
                if (parms$numStrains>1){
                    Bval=sum(stateVarValues[paste0(parms$myPars$mucin.degrader,'.',1:parms$numStrains)])
                }else{
                    Bval=sum(stateVarValues[parms$myPars$mucin.degrader])
                }
                
            }else{

                Bval=0
                for (md in parms$myPars$mucin.degrader){
                    stn=paste0(parms$myPars$mucin.degrader,'.',1:parms$numStrains[md])
                    Bval=Bval+sum(stateVarValues[stn])
                }
                    
                
            }

            lim=Bval/(Bval+Wvol*parms$myPars$halfSat.mucin) #note Kb is in g/L

            
            if (varName=='C'){
                v.host=parms$myPars[['C.from.host']]*lim
            }else if (varName=='RS'){
                v.host=parms$myPars$RS.frac.mucin*parms$myPars[['C.from.host']]*lim
            }else if (varName=='NSP'){
                v.host=(1-parms$myPars$RS.frac.mucin)*parms$myPars[['C.from.host']]*lim
            }else if (varName=='Protein'){
                v.host=parms$myPars[['P.from.host']]*lim
            }
        }else{
            v.host=0
        }

        v=v+v.host
        
    }

    return(v)
}
