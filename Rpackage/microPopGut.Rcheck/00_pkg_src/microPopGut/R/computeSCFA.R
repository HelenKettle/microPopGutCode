#' computeSCFA
#'
#' calculates SCFA (mM) from model results 
#'
#' @param out is the output of microPopGut()$solution[\[n]]. It is a matrix of all model state variables (columns) at each output time (rows) for compartment n 
#' @param SCFAnames vector of strings which are the names of the SCFA
#' @param molarMass is the molar mass of each resource
#' @param waterName is the name of the water variable
#' @export
#' @return returns a matrix with SCFA in each column and rows over time


computeSCFA=function(out,SCFAnames,molarMass,waterName='H2O'){

    L=nrow(out)

    Snames=NULL
    for (s in SCFAnames){
        if (s%in%colnames(out)){
            Snames=c(Snames,s)
        }
    }
    
    TSCFA=rep(NA,L)
    SCFA.mM=matrix(NA,nrow=L,ncol=length(Snames))
    colnames(SCFA.mM)=Snames

    Wvol=out[,waterName]/1000
    
    for (i in 1:nrow(out)){
        for (s in Snames){
            SCFA.mM[i,s]=1000*out[i,s]/molarMass[s]/Wvol[i]
        }
        TSCFA[i]=sum(SCFA.mM[i,Snames])
    }


    return(SCFA.mM)
}
    
