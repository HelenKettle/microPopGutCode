parms=list()
pars=list()
myPars=list()

parms$myPars$waterName='water'

time=1

#growth limit - check it's between 0 and 1
resourceValues=c(R1=10,R2=5)
stateVarValues=c(water=100,R1=10,R2=5)
strainHalfSat=c(R1=0.01,R2=0.01)
strainName=NULL
groupName=NULL
pathName=NULL
varName='R1'
allSubType=c(R1='S',R2='S')
test_that("growth Lim func", {
    expect_lt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R1', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,1)
})
allSubType=c(R1='S',R2='Se')
test_that("growth Lim func", {
    expect_gt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,0)
})
test_that("growth Lim func", {
    expect_lt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,1)
})
allSubType=c(R1='S',R2='Sb')
test_that("growth Lim func", {
    expect_gt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,0)
})
test_that("growth Lim func", {
    expect_lt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,1)
})
allSubType=c(R1='S',R2='Sm')
test_that("growth Lim func", {
    expect_gt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,0)
})
test_that("growth Lim func", {
    expect_lt(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R2', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,1)
})
resourceValues=c(R1=0,R2=5)
stateVarValues=c(water=100,R1=0,R2=5)
test_that("growth Lim func", {
    expect_equal(newRateFuncs$growthLimFunc(strainName, groupName, pathName,
    'R1', resourceValues, allSubType, strainHalfSat, stateVarValues,parms) ,0)
})
