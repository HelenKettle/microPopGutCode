#' newRateFuncs for microPopGut
#'
#' List of functions that are used by the ODE solver
#' these functions can be changed by the user but all must be listed.
#' 
#' rateFuncsDefault=list(pHFunc=pHFuncDefault, pHLimFunc=pHLimFuncDefault, extraGrowthLimFunc=extraGrowthLimFuncDefault, growthLimFunc=growthLimFuncDefault, combineGrowthLimFunc=combineGrowthLimFuncDefault, uptakeFunc=uptakeFuncDefault, productionFunc=productionFuncDefault, combinePathsFunc=combinePathsFuncDefault, massBalanceFunc=massBalanceFuncDefault, entryRateFunc=entryRateFuncDefault, removalRateFunc=removalRateFuncDefault)
#'
#' note that in these functions, the parms list is intrinsic to microPop whereas the myPars list is defined by the inputs to microPopGut() and then added to parms
#' 
#'
#' @include pHFuncMPG.R
#' @include pHLimFuncMPG.R
#' @include entryRateFuncMPG.R
#' @include growthLimFuncMPG.R
#' @include removalRateFuncMPG.R
#'
#' @export
#' 
#'

newRateFuncs=microPop::rateFuncsDefault

newRateFuncs$pHFunc = pHFuncMPG

newRateFuncs$pHLimFunc = pHLimFuncMPG

newRateFuncs$entryRateFunc = entryRateFuncMPG

newRateFuncs$removalRateFunc = removalRateFuncMPG

newRateFuncs$growthLimFunc = growthLimFuncMPG

