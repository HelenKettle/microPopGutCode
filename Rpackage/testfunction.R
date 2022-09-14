mp.orig=function(a,b,c){

    print(c(a,b,c))
    return(c(a,b,c))
}

mp.gut=function(d,e,...){

    a=1
    b=2
    c=3
    
    mpo.args=formalArgs(mp.orig)
    dot.args <- list(...)
    
    for (da in names(dot.args)){
        if (da%in%mpo.args){
            #need to overwrite default arg value
            assign(da,dot.args[da])
        }
    }

    mp.orig(a,b,c)
    
    return(c(d,e))
}

d='D'
e='E'
f='a'
k=mp.gut(d,e,c=7)

#check stuff about lists

pars=list(a=1,b=2,c=3)
meals=list(st=12,fin=14)
new=c(pars,meals)

