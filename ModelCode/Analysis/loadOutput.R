loadOutput=function(fileName,folderName,waterName,SCFAnames=c('A','Bu','Pr'),substrateNames=c('Protein','RS','NSP'),numStrains=NA,fileName2=NA){

#    fileName='numDays10seed1fix.pHFALSEfluc.inflowFALSEfluc.subst.compFALSEnumStrains1BMpd0C.from.diet50P.from.diet10W.from.diet1100transitTime1.25runNumber1RS.frac0.78LPc8.7.Rdata'

#    folderName='/home/helen/MicroPop/Applications/Human/microPopHumanGut/Nov2020/Results/TenGroup/LPsims/'

 #   numStrains=NA
  #  fileName2=NA
#SCFAnames=c('A','Bu','Pr')



    
#     source('computeSCFA.R')

#    print('in loadOutput()')


    
    fn=paste0(folderName,fileName)
        
    if (file.exists(fn)){

        load(fn)

        outs1=allres.list$solution
        nbox=length(outs1)
        
        pars=allres.list$params
        
        ph1=allres.list$pH.vals
        ph.time1=allres.list$pH.time
        
        microbeNames=pars[[1]]$microbeNames
        bioNames=microbeNames

        if (!is.na(numStrains)){
            strainNames=pars[[1]]$allStrainNames
            bioNames=strainNames
        }else{
            strainNames=NULL
        }


        time1=outs1[[1]][,'time']
            
        if (!is.na(fileName2)){
                
            fn2=paste0(folderName,fileName2)
            print(fn2)
            if (file.exists(fn2)){print('there are two parts!')}
            
            load(fn2)
            outs2=allres.list$solution
            ph2=allres.list$pH.vals
            ph.time2=allres.list$pH.time
            time2=outs2[[1]][,'time']+max(time1,na.rm=TRUE)+mean(diff(time1,na.rm=TRUE),na.rm=TRUE)
            time=c(time1,time2)
            
            outs=list()
            ph=list()
            ph.time=list()
            for (box in 1:3){
                outs[[box]]=rbind(outs1[[box]],outs2[[box]])
                ph[[box]]=c(ph1[[box]],ph2[[box]])
                ph.time[[box]]=c(ph.time1[[box]],ph.time2[[box]])
            }
            
        }else{
            time=time1
            ph=ph1
            ph.time=ph.time1
            outs=outs1
        }

        pH=list()
        for (b in 1:nbox){
            if (sum(!is.na(ph[[b]]))>3){
                pH[[b]]=approx(ph.time[[b]],ph[[b]],time)$y 
            }else{
                pH[[b]]=rep(NA,length(time))
            }
        }
  
        
        Lt=length(time)
        microbeNames=pars[[1]]$microbeNames


        SCFAmass.array=array(NA,dim=c(Lt,(length(SCFAnames)+1),nbox),dimnames=list(NULL,c(SCFAnames,'TSCFA'),NULL))
        
        SCFA.array=array(NA,dim=c(Lt,(length(SCFAnames)+1),nbox),dimnames=list(NULL,c(SCFAnames,'TSCFA'),NULL))
        waterMat=matrix(NA,nrow=Lt,ncol=nbox)
        
        for (b in 1:nbox){
            ss=computeSCFA(outs[[b]],SCFAnames,pars[[1]]$molarMass,waterName)
            #print(dim(ss$SCFA))
            SCFA.array[1:nrow(ss$SCFA),1:3,b]=ss$SCFA
            SCFA.array[1:nrow(ss$SCFA),4,b]=ss$TSCFA
            SCFAmass.array[1:nrow(ss$SCFA),1:3,b]=outs[[b]][,SCFAnames]
            SCFAmass.array[1:nrow(ss$SCFA),4,b]=rowSums(outs[[b]][,SCFAnames])
            waterMat[1:nrow(outs[[b]]),b]=outs[[b]][,waterName]
        }

        
        microbe.array=array(NA,dim=c(Lt,length(microbeNames),nbox),dimnames=list(NULL,microbeNames,NULL))
        for (b in 1:nbox){
            mm=outs[[b]][,microbeNames]
            microbe.array[1:nrow(mm),,b]=mm
        }

        substrate.array=array(NA,dim=c(Lt,length(substrateNames),nbox),dimnames=list(NULL,substrateNames,NULL))

        for (b in 1:nbox){
            mm=outs[[b]][,substrateNames]
            substrate.array[1:nrow(mm),,b]=mm
        }


        pH.mat=matrix(NA,nrow=Lt,ncol=nbox)
        for (b in 1:nbox){
            pH.mat[1:length(pH[[b]]),b]=pH[[b]]
        }

        
        final=list(SCFA.array=SCFA.array,microbe.array=microbe.array,ph.mat=pH.mat,
            time=time,pars=pars,water=waterMat,
            SCFAmass.array=SCFAmass.array,substrate.array=substrate.array)

    }else{
        print(paste(fn,'does not exist'))
        final=NA
    }
    
    return(final)

}

        
        
