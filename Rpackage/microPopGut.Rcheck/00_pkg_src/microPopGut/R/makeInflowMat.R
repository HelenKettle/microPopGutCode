makeInflowMat=function (out){
    
    soln = out$solution
    washOutRate = out$parms$Smats$washOut

    mat = soln * NA
    mat[, 1] = soln[, "time"]

    for (i in 2:ncol(soln)) {
        xname = colnames(soln)[i]
        x = microPop::getGroupName(xname, out$parms$microbeNames)
         for (t in 1:nrow(soln)) {
            mat[t, i] = soln[t,xname]*washOutRate[x]
        }
    }
    return(mat)
}
