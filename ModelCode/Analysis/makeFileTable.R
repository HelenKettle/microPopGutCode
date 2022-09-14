makeFileTable=function(files,batch.par.names){
    
    f.table=matrix(NA,ncol=length(batch.par.names),nrow=length(files))
    colnames(f.table)=batch.par.names

    for (f in 1:length(files)){
    
        fn=files[f]

        for (varName in batch.par.names){

    #find start of string
            st=gregexpr(pattern=varName,fn)

            if (st[[1]]>0){
                index.end=unlist(st)+nchar(varName)
                
                varVal=NULL
                for (i in 0:10){
                    
                    s=substr(fn,(index.end+i),(index.end+i))
                    
                    if (s%in%c('.',paste(0:9))){
                        varVal=paste0(varVal,s)
                    }else{
                        break
                    }
                
                }
                
        #test for TRUE or FALSE
                s1=substr(fn,(index.end),(index.end+3))
                if (s1=='TRUE'){varVal=s1}
                
                s2=substr(fn,(index.end),(index.end+4))
                if (s2=='FALSE'){varVal=s2}
                
            #print(paste(varName,'=',varVal))
                f.table[f,varName]=varVal
            }
        }
    }
    return(f.table)
}

