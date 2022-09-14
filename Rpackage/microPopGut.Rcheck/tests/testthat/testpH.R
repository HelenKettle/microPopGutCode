parms=list()
pars=list()
myPars=list()


parms$myPars$SCFAnames=c('A','Bu','Pr')
parms$myPars$waterName='water'
parms$myPars$fix.pH=TRUE
parms$myPars$save.pH=FALSE

parms$molarMass=c('A'=60,'Bu'=74,'Pr'=88)


time=1

stateVarValues=c('water'=10,'A'=10,'Bu'=20,'Pr'=30)

time.ode<<-0

pH.ode<<-6


#--- test fixed pH
parms$myPars$fix.pH=TRUE
parms$myPars$pH.fixed=c(5,6,7)
             
    
for (box in 1:3){

    parms$myPars[['compartment']]=box

    ans=c(5,6,7)[box]
    
    test_that("test 1: pH func (fixed pH)", {
        expect_equal(newRateFuncs$pHFunc(time,parms,stateVarValues),ans)
    })
}


parms$myPars$fix.pH=FALSE

#test pH is 5 for high TSCFA
stateVarValues=c('water'=100,'A'=200,'Bu'=200,'Pr'=200)
test_that("pH func not fixed pH (high TSCFA)", {
    expect_equal(newRateFuncs$pHFunc(time,parms,stateVarValues), 5)
})

stateVarValues=c('water'=100,'A'=0,'Bu'=0,'Pr'=0)
test_that("pH func not fixed pH (low TSCFA)", {
    expect_equal(newRateFuncs$pHFunc(time,parms,stateVarValues), 8)
})

stateVarValues=c('water'=100,'A'=0.3,'Bu'=0.2,'Pr'=0.2)
test_that("pH func not fixed pH (mid TSCFA)", {
    expect_lt(newRateFuncs$pHFunc(time,parms,stateVarValues), 8)
})
test_that("pH func not fixed pH (mid TSCFA)", {
    expect_gt(newRateFuncs$pHFunc(time,parms,stateVarValues), 5)
})

#pH limit func -----

parms$strainPHcorners=rbind(bug1=c(1,2,3,4),bug2=c(3,4,5,6))
pH=2.5
test_that("pH Lim func", {
    expect_equal(newRateFuncs$pHLimFunc('bug1','groupName',2.5,parms),1)
})
test_that("pH Lim func", {
    expect_equal(newRateFuncs$pHLimFunc('bug2','groupName',7,parms),0)
})

