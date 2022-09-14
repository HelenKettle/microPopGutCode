#' pH Function for microPopGut
#'
#' Return the value of pH in pH units
#'
#' @param time (scalar). The current time point in the ODE solver.
#' @param parms List which contains all information required by the ODE solver
#' @param stateVarValues State vector (resources and microbes) (with names)
#' @return (scalar) pH at the given time
#' @export


pHFuncMPG = function (time, parms,stateVarValues){

    
     if (parms$myPars$fix.pH){

         pH=parms$myPars$pH.fixed[ parms$myPars[['compartment']] ]

    }else{

        if (parms$myPars$waterName%in%names(stateVarValues)){
            Wvol=stateVarValues[parms$myPars$waterName]/1000
        }else{
            stop('water is not a state variable so concentrations can not be calculated')
        }
        
        TSCFA=1000*sum(stateVarValues[parms$myPars$SCFAnames]/parms$molarMass[parms$myPars$SCFAnames],na.rm=TRUE)/(Wvol+1e-10)

        pH=8.02-0.0174*TSCFA

        if (pH<5){pH=5}
        if (pH>8){pH=8}
        

 #       if (time>(0.1/60/24)){#after 0.1 minutes...
 #           time.diff=time-max(time.ode)
 #           if (time.diff>(10/24/60)){#save pH every 10 minutes
 #               assign('pH.ode',c(pH.ode,pH),envir=.GlobalEnv)
 #               assign('time.ode',c(time.ode,time),envir=.GlobalEnv)
 #           }
 #       }
    
    }

    pH=unname(pH)
    return(pH)
}
