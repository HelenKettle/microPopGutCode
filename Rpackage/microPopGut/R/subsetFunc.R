subsetFunc = function(array.in, gname, path) {
    
    if (length(dim(array.in)) == 3) {
        v = array.in[gname, , path]
        resnames = colnames(array.in[gname, , path, drop = FALSE])
        names(v) = resnames
    }
    
    if (length(dim(array.in)) == 2) {
        # already subsetted by gname or strain name
        v = as.vector(array.in[path, ])
        resnames = colnames(array.in)
        names(v) = resnames
    }
    
    return(v)
}
