parms=list()
pars=list()
myPars=list()


#entry rate func - note output is a named vector
parms$myPars$compartment=1
varValue=NULL
parms$microbeNames=c('G1','G2')
parms$numStrains=1
inflowRate=c('G1'=10,'G2'=5,'R1'=3,'R2'=6)
parms$myPars$fluc.inflow=FALSE
parms$myPars$substrateNames=c('R1','R2')
parms$myPars$carbs=c('RS','NSP')
stateVarValues=NULL
time=NULL
pars$halfSat.mucin=0.1

test_that("entry ratefunc 1", {
    expect_equal(newRateFuncs$entryRateFunc('G1',varValue,stateVarValues,
                                   time,inflowRate,parms),c(G1=10))
})
test_that("entry ratefunc 2", {
    expect_equal(newRateFuncs$entryRateFunc('G2',varValue,stateVarValues,
                                   time,inflowRate,parms),c(G2=5))
})
#test inflow mat
parms$myPars$compartment=2
parms$myPars$inflow.mat=cbind(time=c(1,2,3),'R1'=c(3,4,5))

time=2
test_that("entry ratefunc 3", {
    expect_equal(newRateFuncs$entryRateFunc('R1',varValue,stateVarValues,
                                   time,inflowRate,parms),c(R1=4))
})

#test mucin
parms$numStrains=1
parms$myPars$compartment=1
parms$myPars$mucin.degrader='G1'
parms$myPars$carbs=c('C','NSP')
inflowRate=c('G1'=10,'G2'=5,'C'=0,'R2'=6)
stateVarValues=c(water=100,C=10,R2=5,G1=10)
parms$myPars$waterName='water'

parms$myPars[['C.from.host']]=10
parms$myPars$halfSat.mucin=0.5

test_that("entry ratefunc 4", {
    expect_gt(newRateFuncs$entryRateFunc('C',varValue,stateVarValues,
                                   time,inflowRate,parms),0)
})

stateVarValues=c(water=100,C=10,R2=5,G1=0,G2=10)
test_that("entry ratefunc 5", {
    expect_equal(newRateFuncs$entryRateFunc('C',varValue,stateVarValues,
                                   time,inflowRate,parms),c(C=0))
})
parms$myPars$mucin.degrader='G2'
test_that("entry ratefunc 6", {
    expect_gt(newRateFuncs$entryRateFunc('C',varValue,stateVarValues,
                                   time,inflowRate,parms),0)
})
