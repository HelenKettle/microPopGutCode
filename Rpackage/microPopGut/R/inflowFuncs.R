#'mealsPattern
#' 
#' this function generates a step function which is 1 for eating time and 0 if not eating. This is then converted to have a mean of 1 over the  model simulation time
#' Note this function is in hours (unlike rest of package)
#' @param num.days number of days simulation is over
#' @param breakfast.start time to start breakfast (24h clock), default is 7am
#' @param lunch.start time to start lunch (24h clock), default is 1pm
#' @param dinner.start time to start dinner (24h clock), default is 7pm
#' @param time.step.h the model time step converted from days to hours
#' @param meal.duration.h length of time for a meal (in hours)
#' @param showPlot Logical produce plot of inflow (defaults to FALSE)

#' @return is a 2 column matrix with time in hours (col 1) and meals (mean=1)
#' @export 
mealsPattern=function(num.days=3,
    breakfast.start=7, lunch.start=13, dinner.start=19,
    time.step.h=1/60, meal.duration.h=0.5, showPlot=FALSE){

    if (lunch.start<breakfast.start){
        stop('lunch.start must be later than breakfast.start. Please use 24 h clock')
    }
    if (dinner.start<breakfast.start | dinner.start<lunch.start){
        stop('dinner.start must be later than breakfast.start and lunch.start. Please use 24 h clock')
    }

    if (meal.duration.h<=time.step.h){
        stop('meal duration must be greater than time step (all in hours)')
    }
    
    times.h=seq(0,num.days*24,by=time.step.h)

    meals=times.h*0

    mod.times.h=times.h%%24

    meals[mod.times.h>=breakfast.start & mod.times.h<(breakfast.start+meal.duration.h)]=1
    meals[mod.times.h>=lunch.start & mod.times.h<(lunch.start+meal.duration.h)]=1
    meals[mod.times.h>=dinner.start & mod.times.h<(dinner.start+meal.duration.h)]=1

                                        #make mean = 1
    if (sum(meals)>0){
        mat=cbind(times.h,meals/mean(meals))
    }else{
        mat=cbind(times.h,meals)
    }

    if (showPlot){
        plot(mat[,1],mat[,2],type='l',xlab='time (h)',ylab='meals')
    }
    
    return(mat)

}

#' mealComposition
    
#' takes step function from mealsPattern() (in hours) then adds in dietary composition for each meal. This is randomised using a gamma distribtion (controlled by scale argument)

#' @param mealsTS output from mealsPattern() 
#' @param substrateNames a vector of substrate names
#' @param substrateMeans a named vector of substrate means, e.g. c(A=1,B=2,C=3)
#' @param scale a named vector controlling the shape of the gamma function, e.g. c(A=0.5,B=1,C=1.5)
#' @param fluc.comp LOGICAL. TRUE for fluctating meal composition
#' @param num.days Scalar Number of days to simulate
#' @param waterName string specifying the name of water e.g. 'H20', 'water' etc
#' @param showPlot FALSE If TRUE will plot meal composition
#' @param saveFig FALSE If TRUE image will be saved as eps
    
    
#' @return output is a matrix with a column for each substrate over time

#' @export
mealComposition=function(mealsTS,substrateNames,substrateMeans,scale,fluc.comp,num.days,waterName,showPlot=FALSE,saveFig=FALSE){


    L=nrow(mealsTS)
    time=mealsTS[,1] #note output from mealsPattern has time in hours
    mag=mealsTS[,2]
    
    Ns=length(substrateNames)
    diet=matrix(0,ncol=Ns,nrow=L)
    colnames(diet)=substrateNames

    
    if (sum(mag)>0){ #i.e. there are some meals occurring in the simulation time

        b=as.numeric(mag>0)
        change.pts=seq(1,L)[diff(b)!=0]

        for (i in seq(1,(length(change.pts)-1),2)){
            f=rep(NA,Ns)
            names(f)=substrateNames
            for (s in substrateNames){
                if (fluc.comp){
                    if (substrateMeans[s]>0){
                        f[s]=stats::rgamma(1,shape=substrateMeans[s]/scale[s],scale=scale[s])
                    }else{
                        f[s]=0
                    }
                }else{
                    f[s]=substrateMeans[s]
                }
            }

            for (j in change.pts[i]:change.pts[i+1]){
                diet[j,]=f
            }
        }
        
        

    if (num.days>=1){
        #rescale so have same daily mean
        for (s in substrateNames){
            if (mean(diet[,s],na.rm=TRUE)>0){
                diet[,s]=substrateMeans[s]*diet[,s]/mean(diet[,s],na.rm=TRUE)
            }else{
                diet[,s]=0*diet[,s]
            }
        }
        
    }else{#if less than a day than set to mean for that fraction of d
        if (max(diet,na.rm=TRUE)>1){
            for (s in substrateNames){
                if (mean(diet[,s],na.rm=TRUE)>0){
                    diet[,s]=num.days*substrateMeans[s]*diet[,s]/mean(diet[,s],na.rm=TRUE)
                }
            }
        }
    }

    
    if (showPlot){
        mealCompBar=matrix(0,ncol=length(substrateNames),nrow=length(change.pts))
        colnames(mealCompBar)=substrateNames
        
        for (s in substrateNames){
            ct=1
            for (j in 1:(length(change.pts)-1)){
                ii=(seq(1,nrow(diet))>change.pts[j] & seq(1,nrow(diet))<change.pts[j+1])
                mealCompBar[ct,s]=mean(diet[ii,s])
                ct=ct+1
            }
        }
                                        #for plotting make into meal blocks
        dev.new(width=10,height=6)
        par(mfrow=c(1,2))
        ii=mealCompBar[,1]>0
        if (sum(ii)>1){
            barplot(t(mealCompBar[ii,substrateNames[1:3]]),xlab=paste('meals over',round(max(time)/24,2),'d'),ylab='substrate (g/d)',col=1:3,ylim=c(0,1.2*max(rowSums(mealCompBar[ii,substrateNames[1:3]]),na.rm=TRUE)))
        
            legend('topright',col=1:3,pch=15,legend=substrateNames[1:3],bty='n')
        }else{
            barplot(t(mealCompBar[ii,substrateNames[1:3]]),xlab=paste('meals over',round(max(time)/24,2),'d'),ylab='substrate (g/d)',col=1:3)
        }
        barplot(t(mealCompBar[ii,waterName]),xlab=paste('meals over',round(max(time)/24,2),'d'),ylab=paste(waterName,'(g/d)'),col='blue')
        
        if (saveFig){
            dev.copy2eps(file='mealComp.eps')
        }
        
    }
    }
    
    return(diet)
}



smallIntestine=function(meals,time.to.reach.colon.h,showPlot=FALSE){

    #' smallIntestine
    #' put output from mealComposition(), combined with a time vec (col 1) through a one compartment ODE model (representing stomach and small intestine)

    #' @param meals output from mealComposition()
    #' @param time.to.reach.colon.h time in hours taken for food to reach the entrance to the colon
    #' @param showPlot Default is FALSE
    
    # @return is a matrix with [time (days), inflow to colon (no units (mean=1))

    #' @importFrom deSolve ode
    #' @importFrom stats approx
    #' @export
    if (time.to.reach.colon.h<=0){
        stop('time.to.reach.colon.h must be greater than 0')
    }
    
    init=c(S=0)

    deriv=function(t,y,parms){
    
        S=y[1]
        S.in=approx(meals[,1],meals[,2],t,rule=2)$y
        dS=S.in-V*S
        
        return(list(dS))
    }

    V=1/time.to.reach.colon.h

    out=deSolve::ode(init,meals[,1],func=deriv,parms=NULL)

    times.d=out[,1]/24

    outflow=V*out[,2]
    
#    if (showPlot){
#        plot(times.d,outflow,type='l',xlab='time (d)',ylab='meals inflow to colon')
#    }
  
    return(cbind(times.d,outflow))
}


#' makeSubstrateMat

#' this function calls mealsPattern(), mealComposition() and smallIntestine()

#' @param substrateNames a vector of substrate names
#' @param substrateMeans a named vector of substrate means, e.g. c(A=1,B=2,C=3)
#' @param num.days number of days simulation is over
#' @param breakfast.start time to start breakfast (24h clock), default is 7am
#' @param lunch.start time to start lunch (24h clock), default is 1pm
#' @param dinner.start time to start dinner (24h clock), default is 7pm
#' @param time.step.h the model time step converted from days to hours
#' @param meal.duration.h length of time for a meal (in hours)
#' @param fluc.comp LOGICAL. TRUE for fluctating meal composition
#' @param time.to.reach.colon.h time from ingestion to reaching colon (hours)
#' @param gamma.mag a number controlling the shape of the gamma function
#' @param RS.frac Fraction of carbs that is RS
#' @param showPlot produce plot of inflow
#' @param saveFig save plot of inflow
#' @param waterName string for the name of the water variable. Default is 'H2O'
#' @param waterMean number for mean inflow of water to colon (g/d). Default is 1100 g/d

#' @return a matrix for substrate inflow (g/d) for each substrate and water (cols) at each model time step (rows). The first col is time in days.
    


#' @import graphics
#' @import grDevices
#'
#' @export
makeSubstrateMat=function(
    substrateNames,
    substrateMeans,
    num.days=2,
    breakfast.start=7,
    lunch.start=13,
    dinner.start=19,
    time.step.h=0.1,
    meal.duration.h=0.5,
    fluc.comp=TRUE,
    time.to.reach.colon.h=7,
    gamma.mag=1,
    RS.frac=0.78,
    showPlot=TRUE,
    saveFig=FALSE,
    waterName='H2O',
    waterMean=1100
){


    #Note, substrateNames input does not include water
    if (gamma.mag<=0){stop('gamma must be greater than zero')}

    allNames=c(substrateNames,waterName)
    allMeans=c(substrateMeans,waterMean)
    names(allMeans)=allNames
    
    gamma=gamma.mag*c(0.5*substrateMeans[substrateNames],0.1*waterMean)
    names(gamma)=allNames
    
    Ns=length(allNames)
    
    #get step function of 0 (no food) and 1 (meal)
    mealsTS=mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,
        time.step.h,meal.duration.h)
    #note mealsTS has time in hours
 
    time=mealsTS[,1] #note output from mealsPattern has time in hours
    mag=mealsTS[,2]
    
    diet=mealComposition(mealsTS,allNames,allMeans,gamma,fluc.comp,num.days,waterName,showPlot,saveFig)

    
    #check numbers are ok before small intestine!
    print('required means are:')
    print(allMeans)
    print(paste('over',num.days,'days, fluctuating means are:'))
    for (sn in allNames){
        print(paste(sn,round(mean(diet[,sn]),2)))
    }

    
    final.diet=diet*0

    #pass through small intestine
    if (sum(mag)>0){
        print('start pass through small intestine')
        for (s in allNames){
            #print(s)
            if (allMeans[s]>0){
                final.diet[,s]=smallIntestine(cbind(mealsTS[,1],diet[,s]),time.to.reach.colon.h)[,2]
            }else{
                final.diet[,s]=0*diet[,s]
            }
        }
        print('finish pass  through small intestine')
    }
    
 
    if (showPlot){
        dev.new(width=10,height=5)
        par(mfrow=c(1,4),mar=c(4,5,2,1))
        for (s in allNames){
            plot(time/24,final.diet[,s],
                 type='l',
                 xlab='time (d)',ylab='inflow rate (g/d)',
                 main=paste(s),cex.lab=1.5)
        }
        if (saveFig){
            dev.copy2eps(file='outflowFromSmallIntestine.eps')
        }
    }
    
    #check numbers are ok!
    print('means are:')
    for (sn in allNames){
        print(paste(sn,round(mean(final.diet[,sn]),2)))
    }

    #output has time in days and a each column is a substrate (and water)
    return(cbind('time'=time/24,final.diet))
    
}
