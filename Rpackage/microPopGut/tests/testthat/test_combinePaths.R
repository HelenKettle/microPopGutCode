
#combining paths
growthRate=c(10,4,6)
test_that("combine paths",{
    expect_equal(newRateFuncs$combinePathsFunc(strainName=NULL,groupName=NULL, growthRate, num.paths=3, pathNames=c('p1','p2','p3')),c('p1'=0.5,'p2'=0.2,'p3'=0.3))
})
growthRate=c(0,-1,0)
test_that("combine paths",{
    expect_equal(newRateFuncs$combinePathsFunc(strainName=NULL,groupName=NULL, growthRate, num.paths=3, pathNames=c('p1','p2','p3')),c('p1'=1,'p2'=01,'p3'=1))
})
