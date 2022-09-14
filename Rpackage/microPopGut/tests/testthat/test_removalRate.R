parms=list()
pars=list()
myPars=list()


parms$myPars$SCFAnames=c('A','Bu','Pr')
parms$myPars$waterName='water'


#removal rate - note output is a named vector
parms$myPars$BMpd=0
varName='R1'
stateVarValues=c(water=100,C=10,R1=5,G1=10)
varValue=stateVarValues[varName]
time=NULL
washOut=c(water=1,C=1,R1=1,G1=1)

parms$myPars$SCFAnames=c('A','Bu','Pr')
test_that("removal rate 1",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),c(R1=5))
})


stateVarValues=c(water=100,C=10,R1=0,G1=10)
varValue=stateVarValues[varName]

test_that("removal rate 2",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),0)
})

stateVarValues=c(water=100,C=10,R1=10,G1=10)
varValue=stateVarValues[varName]
washOut=c(water=1,C=1,R1=10,G1=1)

test_that("removal rate 3",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),c(R1=100))
})


#absorption
parms$myPars$abs.scfa=10
parms$myPars$SCFAnames=c('A','Bu','Pr','R1')
washOut=c(water=1,C=1,R1=0,G1=1)
test_that("removal rate 4",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),c(R1=100))
})
washOut=c(water=1,C=1,R1=10,G1=1)
test_that("removal rate",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),c(R1=200))
})

parms$myPars$abs.water=10
varName='water'
stateVarValues=c(water=100,C=10,R1=5,G1=10)
varValue=stateVarValues[varName]
washOut=c(water=0,C=1,R1=10,G1=1)

test_that("removal rate",{
    expect_equal(newRateFuncs$removalRateFunc(varName, varValue, stateVarValues, time, washOut, parms),c(water=1000))
})
