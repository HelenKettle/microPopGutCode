#' microPopGut

#' runs the human colon version of microPop

#' @param numDays number of days over which to run simulation
#' @param time.step simulation time step (in days)
#' @param transitTime time taken to travel through colon (d)
#' @param microbeNames named vector with names of microbial data frames
#' @param microbeNames.short named vector with shortened names for plotting 
#' @param numStrains scalar if all groups have the same number of strains or a named vector specifying the number of strains in each microbial group. 

#' @param absorption named vector with specific absorption rates (/d) for water and SCFA. Default is c('water'=3,'SCFA'=9.6)

#' @param pH.pars list controlling pH in the model. Default is list(pH.fixed=c(5.7,6.2,6.6),fix.pH=FALSE,pHLimit=TRUE),
#' \itemize{
#' \item    pH.fixed=c(5.7,6.2,6.6) pH values in prox, trans and distal compartments - used if fix.pH is TRUE
#' \item    fix.pH=FALSE Set to TRUE to keep pH fixed. If FALSE, pH varies with SCFA concentrations
#' \item    pHLimit=TRUE If TRUE then microbial growth is pH limited
#'}

#' @param colon.pars list of parameters describing colon. Default is list(vol.litres=c(0.41,0.98,1.63),compartmentNames=c('prox','trans','dist'))

#' @param bowel.movements List of parameters for modelling bowel movements. Default is list(BM.duration.h=15/60,frac.distal.emptied=0.90, BMpd=0,start.BM.time=list(7,c(7,19),c(7,15,21))
#' \itemize{
#' \item BM.duration.h time for bowel movement in hours
#' \item frac.distal.emptied fraction of the distal colon that is emptied by bowel movemen
#' \item BMpd number of bowel movements per day
#' \item start.BM.time list of times of bowel movements for 1, 2, or 3 BM per day
#' }

#' @param init list of initial values for carbohydrates (C), protein (P), microbes (B), Acetate, Propionate and Butyrate and water (W) in grams. Default values are: list(C=1,P=1,B=10,Acetate=0.3606,Propionate=0.1482,Butyrate=0.1762, W=100). SCFA initial values are based on concentrations of 60 mM, 20 mM, 20 mM for acetate, propionate and butyrate. The microbial mass is split evenly between the microbial groups. The initial conditions are same in each model compartment.

#' @param inflow list specifying the inflow of substrates (C, P) and water (W) in g/d, to the colon. Default values are: list(C=50,P=10,W=1100,RS.frac=0.78). The RS.frac is the fraction of RS in the carbohydrate i.e. RS/(RS+NSP) assuming RS+NSP=C
    
#' @param meals list of parameters describing meals. Default is list(seed=1,fluc.inflow=FALSE,fluc.subst.comp=FALSE,plotInflow=FALSE,saveInflowFig=FALSE,breakfast.start=7,lunch.start=13,dinner.start=19,meal.duration.h=0.5,time.to.reach.colon.h=7,gamma.mag=1)
#' \itemize{
#' \item seed seed for random number generator used for stochastic meal composition only (fluc.subst.comp=TRUE)
#' \item fluc.inflow FALSE if incoming substrate is constant i.e. TRUE means meals
#' \item fluc.subst.comp FALSE change to TRUE for substrate composition varying randomly
#' \item plotInflow FALSE change to TRUE to plot the incoming substrate
#' \item saveInflowFig FALSE Change to TRUE to save inflow figs
#' \item breakfast.start time breakfast is eaten (24 h clock)
#' \item lunch.start time lunch is eaten (24 h clock)
#' \item dinner.start time dinner is eaten (24 h clock)
#' \item meal.duration.h time taken to eat a meal (in hours)
#' \item time.to.reach.colon.h time taken to move through stomach and small intestine etc (h)
#' \item gamma.mag (positive number) this determines the level of variation in diet fluctuations (1.1 will give bigger fluctuations, 0.9 smaller).
#' }

#' @param
#' mucin list of parameters describing mucins. Default is list(mucin.rate=5,RS.frac.mucin=0.5,mucin.degrader=c('Bacteroides','NoButyFibreDeg'),halfSat.mucin=0.5,P.frac.mucin=0.2,C.frac.mucin=0.75)
#' \itemize{
#' \item mucin.rate  rate of mucin (g/d) secreted into colon
#' \item RS.frac.mucin fraction of C in mucin that is RS
#' \item mucin.degrader microbes that break down mucin
#' \item halfSat.mucin half saturation constant (g/l of B) for mucin breakdown
#' \item P.frac.mucin fraction of mucin that is protein
#' \item C.frac.mucin fraction of mucin that is carbohydrate
#'}

#' @param odeFunc The function the ODE solver in microPopModel() will use - the default is derivsDefault provided by the package but if the user wants to make significant changes a new ODE function file can be used. See ?derivsDefault
#' @param oneStrainRandomParams Logical to allow randomization of parameters even if there is only one strain. The default is FALSE which means that if numStrains=1 then the group params are used; if numStrains>1 then the parameters are automatically randomised according to info given in strainOptions. If oneStrainRandomParams=TRUE then even if there is only one strain its parameters will be randomised according to info given in strainOptions. 
#' @param plotOptions List containing instructions for plotting: Default is list(plotFig=TRUE, sumOverStrains=FALSE, saveFig=FALSE, figType='eps', figName='microPopFig', yLabel='Concentration (g/L)', xLabel='Time').\cr
#' To turn off plot generation set plotFig=FALSE. If there are multiple strains these are all plotted if sumOverStrains=FALSE, otherwise they will be summed over each group. To save plot, saveFig=TRUE, figType (format) can be 'eps','png', 'pdf' or 'tiff' and is specified in figType (string), the name is figName (string) to which the string 'Microbes' or 'Resources' will be added for the respective plots.
#' @param odeOptions List containing instructions for the ODE solver ('deSolve'). Default: list('atol'=1e-6,'rtol'=1e-6,'method'='lsoda'). See ?ode for more details.
#' @param strainOptions List containing instructions for specifying strain parameters. Default: list(randomParams=c('halfSat', 'yield', 'maxGrowthRate', 'pHtrait'), seed=1, distribution='uniform', percentTraitRange=0, maxPHshift=0., applyTradeOffs=FALSE, tradeOffParams=NULL, paramsSpecified=FALSE, paramDataName=NULL). 
#' \itemize{
#' \item randomParams (vector) specifying which parameters need to be stochastically generated.
#' \item seed (number) seed for random number generator.
#' \item distribution (string) - either 'uniform' or 'normal' specifying the shape of the distribution from which to draw the random strain parameters. 
#'\item percentTraitRange (number) this is the percentage either side of the group parameter value which the strain parameter may range e.g. if percentTraitRange=10 then range is 0.9x to 1.1x for group mean x. 
#'\item maxPHshift (number) pH units to range over.
#' \item applyTradeOffs (logical) to trade off `good' and `bad' parameter values. 
#'\item tradeOffParams (vector of two strings) - parameters to trade off against each other. Note that pHtrait can not be traded off as whether this trait is good or bad depends on the environmental pH.
#'\item paramsSpecified (logical) TRUE if strain parameters are read in from a file (whose name is specified in paramDataName). The file must have colnames c(strainName, paramName, paramVal, paramUnit, resource,path) and where strainName is in format 'groupName.i' where i is the strain number.
#' }
#' @param checkingOptions (List) Default is list(checkMassConv=FALSE, balanceTol=1e-2, reBalanceStoichiom=FALSE, stoiTol=0.1, checkForNegs=TRUE, negTol=-1e-2). 
#' \itemize{
#' \item checkMassConv=TRUE checks for mass conservation in the ODE solver with a tolerance of 'balanceTol' (default is FALSE).
#' \item reBalanceStoichiom will check the mass balance of the stoichiometries on every metabolic path and rebalance if these are not conserving mass within a tolerance of stoiTol (a warning message will be issued). Rebalancing will only affect the final solution if the pathway contains only essential resources (Rtype 'Se') and microbial biomass is a product (Rtype 'Pb').  
#' \item checkForNegs If TRUE the function checkSolution is called and the solution for each variable, x, is checked for negative values that are greater in magnitude than negTol*max(x). If negative values occur then the solution is incorect and either the problem is incorrectly specified or the tolerances in the ODE solver need to be smaller.
#' }
#' @param microbeMolarMass Scalar. Mass of 1 mole of microbes - default is 113g/mol (Batstone et al., 2002)
#' @param bacCutOff Scalar. Amount of bacteria below which the bacteria are considered to have left the system and can't grow, default =1e-14. If this is set to zero then bacteria will always be able to grow again as zero is never reached.
#' @param ... Add your own inputs

#' @return The output is a list containing lists called 'solution', 'params' and myPars. These lists have 3 components (one for each model compartment (i.e. the proximal, transverse and distal colon). Solution[\[n]] (for n=1,2,or 3 depending on compartment) is a matrix where rows are points in time and the columns are the state variables. params[\[n]] contains all the information needed to run the model. Use names(out$params[\[1]]) to see contents.

#' @import microPop

#' @export
microPopGut=function(
    numDays,
    time.step,
    transitTime,
    microbeNames,
    microbeNames.short,
    numStrains=1,
    
    absorption=c(
        'water'=3,
        'SCFA'=9.6),

    pH.pars=list(
        pH.fixed=c(5.7,6.2,6.6),
        fix.pH=FALSE,
        pHLimit=TRUE),

    colon.pars=list(
        vol.litres=c(0.41,0.98,1.63),
        compartmentNames=c('prox','trans','dist')),

    bowel.movements=list(
        BM.duration.h=15/60,
        frac.distal.emptied=0.90,
        BMpd=0, 
        start.BM.time=list(7,c(7,19),c(7,15,21))), 

    init=list(
        C=1,
        P=1,
        B=10,
        Acetate=0.3606,
        Propionate=0.1482,
        Butyrate=0.1762,
        W=100),

    inflow=list(
        C=50,
        P=10,
        W=1100, 
        RS.frac=0.78),     

    meals=list(
        seed=1, 
        fluc.inflow=FALSE,
        fluc.subst.comp=FALSE,
        plotInflow=FALSE,
        saveInflowFig=FALSE,
        breakfast.start=7,
        lunch.start=13, 
        dinner.start=19,
        meal.duration.h=0.5, 
        time.to.reach.colon.h=7,
        gamma.mag=1 
        ),
    
    mucin=list(mucin.rate=5,
        RS.frac.mucin=0.5,
        mucin.degrader=c('Bacteroides','NoButyFibreDeg'),
        halfSat.mucin=0.5, 
        P.frac.mucin=0.2, 
        C.frac.mucin=0.75
    ),
    
    #defaults to pass in microPopModel()
    odeFunc=microPop::derivsDefault,
    oneStrainRandomParams = FALSE,
    microbeMolarMass = 113,
    bacCutOff = 1e-14,

    plotOptions = list(yLabel = "Concentration (g/L)", xLabel = "Time",
        plotFig = FALSE, 
        sumOverStrains = FALSE, saveFig = FALSE, figType = "eps", figName = "microPopFig"),
    
    odeOptions = list(atol = 1e-06, rtol = 1e-06, method = "lsoda"),
    
    strainOptions = list(randomParams = c("halfSat", "yield", "maxGrowthRate", 
        "pHtrait"), seed = 1, distribution = "uniform", percentTraitRange = 0,
        maxPHshift = 0, 
        applyTradeOffs = FALSE, tradeOffParams = NULL, paramsSpecified = FALSE,
        paramDataName = NULL),
    
    checkingOptions = list(checkMassConv = FALSE, balanceTol = 0.01,
        reBalanceStoichiom = FALSE, 
        stoiTol = 0.1, checkForNegs = TRUE, negTol = -0.01, checkStoichiomBalance = TRUE),
  
    ...){

#----------------------------------------------------------
    absorption.default=c('water'=3,'SCFA'=9.6)

    pH.pars.default=list(pH.fixed=c(5.7,6.2,6.6),fix.pH=FALSE, pHLimit=TRUE)

    colon.pars.default=list(vol.litres=c(0.41,0.98,1.63),
        compartmentNames=c('prox','trans','dist'))

    bowel.movements.default=list(BM.duration.h=15/60, frac.distal.emptied=0.90,
        BMpd=0,  start.BM.time=list(7,c(7,19),c(7,15,21)))

    init.default=list(C=1, P=1, B=10, Acetate=0.3606, Propionate=0.1482,
        Butyrate=0.1762,  W=100)

    inflow.default=list( C=50, P=10,  W=1100, RS.frac=0.78)     

    meals.default=list(seed=1, fluc.inflow=FALSE, fluc.subst.comp=FALSE,
        plotInflow=FALSE,saveInflowFig=FALSE,
        breakfast.start=7, lunch.start=13,  dinner.start=19,
        meal.duration.h=0.5, 
        time.to.reach.colon.h=7,  gamma.mag=1 )

    mucin.default=list(mucin.rate=5, RS.frac.mucin=0.5,
        mucin.degrader=c('Bacteroides','NoButyFibreDeg'),
        halfSat.mucin=0.5, P.frac.mucin=0.2, C.frac.mucin=0.75)

    
    absorption = replaceListItems(absorption, absorption.default)
    pH.pars. = replaceListItems(pH.pars,pH.pars.default)
    colon.pars = replaceListItems(colon.pars,colon.pars.default)
    bowel.movements = replaceListItems(bowel.movements,bowel.movements.default)
    init = replaceListItems(init,init.default)
    inflow = replaceListItems(inflow,inflow.default)
    meals = replaceListItems(meals,meals.default)
    mucin = replaceListItems(mucin,mucin.default)

#---------------------------------------------------------------------------    
    pars=list()
    
    pars$numDays=numDays
    pars$time.step=time.step
    pars$transitTime=transitTime
    pars$numStrains=numStrains
    pars$microbeNames=microbeNames
    pars$microbeNames.short=microbeNames.short
    
    #absorption
    pars$abs.scfa=absorption['SCFA'] #9.6, #/d Ruppin et al 1980
    pars$abs.water=absorption['water'] 

    #colon
    pars$prox.frac=colon.pars$vol.litres[1]/sum(colon.pars$vol.litres)
    pars$trans.frac=colon.pars$vol.litres[2]/sum(colon.pars$vol.litres)
    pars$dist.frac=colon.pars$vol.litres[3]/sum(colon.pars$vol.litres)
    pars$compartmentNames=colon.pars$compartmentNames

    #init
    for (nm in names(init)){
        nnm=paste0(nm,'0')
        pars[[nnm]]=init[[nm]]
    }

        #inflow
    for (nm in names(inflow)){
        if (nm!='RS.frac'){
            nnm=paste0(nm,'.from.diet')
        }else{
            nnm=nm
        }
        pars[[nnm]]=inflow[[nm]]
    }          

    pars=c(pars,pH.pars,bowel.movements,inflow,meals,mucin)
    
    #defaults to pass in microPopModel()
    plotOptions.default = list(yLabel = "Concentration (g/L)", xLabel = "Time",
        plotFig = FALSE, 
        sumOverStrains = FALSE, resourceLegendPosition="topleft",
        microbeLegendPosition="topleft",
        saveFig = FALSE, figType = "eps", figName = "microPopFig")
    
    odeOptions.default = list(atol = 1e-06, rtol = 1e-06, method = "lsoda")
    
    strainOptions.default = list(randomParams = c("halfSat", "yield", "maxGrowthRate", 
        "pHtrait"), seed = 1, distribution = "uniform", percentTraitRange = 0,
        maxPHshift = 0, 
        applyTradeOffs = FALSE, tradeOffParams = NULL, paramsSpecified = FALSE,
        paramDataName = NULL)
    
    checkingOptions.default = list(checkMassConv = FALSE, balanceTol = 0.01,
        reBalanceStoichiom = FALSE, 
        stoiTol = 0.1, checkForNegs = TRUE, negTol = -0.01,
        checkStoichiomBalance = TRUE)
    

    
    plotOptions = replaceListItems(plotOptions, plotOptions.default)
    odeOptions = replaceListItems(odeOptions, odeOptions.default)
    strainOptions = replaceListItems(strainOptions, strainOptions.default)
    checkingOptions = replaceListItems(checkingOptions, checkingOptions.default)
 
    
    #add in fixed parameters
    pars$nbox=3 #number of model compartments 
    pars$SCFAnames=c('Acetate','Butyrate','Propionate')
    pars$carbs=c('RS','NSP')
    pars$waterName='H2O'
    resSI=microPop::resourceSysInfoHuman #intrinsic to microPop
    micSI=microPop::microbeSysInfoHuman #intrinsic to microPop
    #---------------------------------

    #checks
    
    #bowel movements
    if (pars$BMpd>0){ 
        pars$V.bm=-log(1-pars$frac.distal.emptied)/(pars$BM.duration.h/24)
        #print(paste('V.BM=',pars$V.bm,'/d'))
        pars$transitTime=(1/pars$BMpd)/pars$dist.frac
        start.BM.time=pars[[pars$BMpd]]
    }

    pars$RS.from.diet=pars$RS.frac*pars$C.from.diet
    pars$NSP.from.diet=(1-pars$RS.frac)*pars$C.from.diet


    substrateMeans=c('Protein'=pars$P.from.diet,'RS'=pars$RS.from.diet,'NSP'=pars$NSP.from.diet)
    substrateNames=names(substrateMeans)
    pars$substrateNames=substrateNames
    waterName=pars[['waterName']]



#INFLOW/OUTFLOW
#note washOut is set via pars$transitTime & and inflow by dietary info depending on compartment
#therefore just set to zero here

    micSI['inflowRate',]=0
    micSI['washOut',]=0
    resSI['inflowRate',]=0
    resSI['washOut',]=0

#INITIAL CONDITIONS (in res and mic dfs) 
    resSI['startValue',c('RS','NSP')]=c(pars$RS.frac,(1-pars$RS.frac))*as.numeric(pars['C0'])
    resSI['startValue',waterName]=pars['W0']
    resSI['startValue','Protein']=pars['P0']
    Wvol0=resSI['startValue',waterName]/1000
    resSI['startValue',pars$SCFAnames]=pars[paste0(pars$SCFAnames,'0')]
    micSI['startValue',pars$microbeNames]=as.numeric(pars['B0'])/length(pars$microbeNames) #start with B0 g evenly spread


    num.time.steps=floor(pars$numDays/pars$time.step)
    simulation.times=seq(0,pars$numDays,length=num.time.steps)

    time.to.reach.colon.h=pars$time.to.reach.colon.h
    
#mucin degraders - check they are present
    m.degraders=NULL
    for (m in pars$mucin.degrader){
        if (m%in%pars$microbeNames){
            m.degraders=c(m.degraders,m)
        }
    }
    pars$mucin.degrader=m.degraders
    
#---------------------------------------------------------------------------

    outs=list()
    p.list=list()
    pH.list=list()
    pH.time.list=list()
    myPars.list=list()

    #assign('resSI',get('resSI'),envir=.GlobalEnv)
    #assign('micSI',get('micSI'),envir=.GlobalEnv)


    for (compartment in 1:pars$nbox){

        print(paste('simulating growth in compartment',compartment,'- please wait!'))

        pars$compartment=compartment
        
        if (compartment==1){
            
            resSI['inflowRate','C']=pars$C.from.diet
            resSI['inflowRate','Protein']=pars$P.from.diet
            resSI['inflowRate',waterName]=pars$W.from.diet
            
            if ('RS'%in%colnames(resSI)){resSI['inflowRate','RS']=pars$RS.from.diet}
            if ('NSP'%in%colnames(resSI)){resSI['inflowRate','NSP']=pars$NSP.from.diet}
            
            resSI['washOut',]=1/(pars$prox.frac*pars$transitTime)
            micSI['washOut',]=1/(pars$prox.frac*pars$transitTime)
                                        #divide mucin (g/d) over length of colon
            mucin.frac=pars$prox.frac
            
            if (pars$fluc.inflow){
                time.step=mean(diff(simulation.times))
                
                subst.mat=makeSubstrateMat(substrateNames,
                    substrateMeans=substrateMeans,
                    num.days=pars$numDays,
                    breakfast.start=pars$breakfast.start,
                    lunch.start=pars$lunch.start,
                    dinner.start=pars$dinner.start,
                    time.step.h=time.step*24,
                    meal.duration.h=pars$meal.duration.h,
                    fluc.comp=pars$fluc.subst.comp,
                    time.to.reach.colon.h=pars$time.to.reach.colon.h,
                    gamma.mag=pars$gamma.mag,
                    RS.frac=pars$RS.frac,
                    showPlot=pars$plotInflow,
                    saveFig=pars$saveInflowFig,
                    waterName=waterName,
                    waterMean=pars$W.from.diet
                                           )
            
                pars$subst.inflow=subst.mat
                
            }
        
        }else if (compartment==2){
            out1=out
            pars$inflow.mat=makeInflowFromSoln(out1)
            resSI['washOut',]=1/(pars$trans.frac*pars$transitTime)
            micSI['washOut',]=1/(pars$trans.frac*pars$transitTime)
            mucin.frac=pars$trans.frac
            
        }else{
            out2=out
            pars$inflow.mat=makeInflowFromSoln(out2)
            resSI['washOut',]=1/(pars$dist.frac*pars$transitTime)
            micSI['washOut',]=1/(pars$dist.frac*pars$transitTime)
            mucin.frac=pars$dist.frac
        }
        
        
        pars$P.from.host=pars$P.frac.mucin*mucin.frac*pars$mucin.rate#g/d
        pars$C.from.host=pars$C.frac.mucin*mucin.frac*pars$mucin.rate #g/d

        
#        assign('pH.ode',NA,envir=.GlobalEnv)
#        assign('time.ode',0,envir=.GlobalEnv)


        out=microPopModel(
            microbeNames=pars$microbeNames,
            times=simulation.times,
            resourceSysInfo=resSI,
            microbeSysInfo=micSI,
            pHLimit=pars$pHLimit,
            numStrains=pars$numStrains,
            rateFuncs=newRateFuncs,
            strainOptions=strainOptions,
            plotOptions=plotOptions,
            checkingOptions=checkingOptions,
            oneStrainRandomParams=oneStrainRandomParams,
            microbeMolarMass=microbeMolarMass,
            bacCutOff=bacCutOff,
            myPars=pars,
            gutModel=TRUE
                                        #odeOptions=list(atol=1e-8,method='lsoda')
        )

#        print('finished microPopModel')

#        print(pH.ode)
 #       solnmat=cbind(out$solution,approx(time.ode,pH.ode,out$solution[,1],rule=2)$y)
        
#        outs[[compartment]]=solmat
        #print(out$solution[1,])
        outs[[compartment]]=out$solution
        p.list[[compartment]]=out$parms
        
        
#        pH.list[[compartment]]=pH.ode#approx(time.ode,pH.ode,simulation.times)$y
#        pH.time.list[[compartment]]=time.ode

        myPars.list[[compartment]]=pars
        
        
    }

    
#    allres.list=list(solution=outs,params=p.list,pH.vals=pH.list,pH.time=pH.time.list,myPars=myPars.list)

    allres.list=list(solution=outs,params=p.list,myPars=myPars.list)

    return(allres.list)

}
