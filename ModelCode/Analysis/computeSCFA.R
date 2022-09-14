computeSCFA=function(out,SCFAnames,molarMass,waterName='H2O'){

    #out is a matrix with rows over time and cols of state variables
    
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


    return(list(SCFA=SCFA.mM,TSCFA=TSCFA))
}
    
