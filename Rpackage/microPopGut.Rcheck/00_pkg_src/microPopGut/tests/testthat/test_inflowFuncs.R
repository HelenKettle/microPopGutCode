#rm(list=ls())
#graphics.off()

#library(testthat)
#source('../../R/inflowFuncs.R')

#in inflowFuncs.R the function makeSubstrateMat() calls mealsPattern(), mealComposition() and smallIntestine()

num.days=2
breakfast.start=7
lunch.start=13
dinner.start=19
time.step.h=1/60
meal.duration.h=0.5

#dev.new()
#par(mfrow=c(2,3),mar=c(5,5,2,1))

expect_error(mealsPattern(num.days,breakfast.start=14,lunch.start,dinner.start,time.step.h,meal.duration.h))
expect_error(mealsPattern(num.days,breakfast.start,lunch.start=20,dinner.start,time.step.h,meal.duration.h))


#mean is 1
test_that("mealsPattern", {
    expect_equal(mean(mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h)[,2]),1)
})

#no values less than zero in either column
test_that("mealsPattern", {
    expect_false(any(mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h)[,2]<0))
})
test_that("mealsPattern", {
    expect_false(any(mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h)[,1]<0))
})

#check is a step function i.e. only 2 unique values
test_that("mealsPattern", {
    expect_length(unique(mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h)[,2]),2)
})
#check one of the unique values is zero
test_that("mealsPattern", {
    expect_true(0%in%unique(mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h)[,2]))
})

#look at plot
meals=mealsPattern(num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,TRUE)

total.food1=sum(meals[,2]*time.step.h)

#mealComposition------------------------------
substrateNames=c('A','B','C')
substrateMeans=c(A=1,B=2,C=3)
scale=c(A=0.5,B=1,C=1.5)

test_that("mealsComposition", {
    expect_equal(colMeans(mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=FALSE,num.days)),substrateMeans)
})
test_that("mealsComposition", {
    expect_equal(colMeans(mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=TRUE,num.days)),substrateMeans)
})

for (s in substrateNames){
    test_that("mealsComposition", {
        expect_length(unique(mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=FALSE,num.days)[,s]),2)
    })
}

for (s in substrateNames){
    test_that("mealsComposition", {
        expect_length(unique(mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=FALSE,num.days)[,s]),2)
    })
    test_that("mealsComposition", {
        expect_true(0%in%unique(mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=FALSE,num.days)[,s]))
    })
}

meals3=mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=FALSE,num.days)

#matplot(meals3,type='l')
meals4=mealComposition(meals,substrateNames,substrateMeans,scale,fluc.comp=TRUE,num.days)

#matplot(meals4,type='l')

#small intestine-------------

time.to.reach.colon.h=7
#look at plot
meals2=smallIntestine(meals,time.to.reach.colon.h,showPlot=FALSE)

#check all values (time and flow) are not less than zero
test_that("smallIntestine", {
    expect_false(any(smallIntestine(meals,time.to.reach.colon.h)<0))
})

total.food2=sum(meals2[,2]*time.step.h)

#expect to have less food than before as all is not out of small intestine yet
test_that('total food',{
    expect_lte(total.food2,total.food1)
})

#make substrate mat------------
gamma.mag=1
RS.frac=0.78
makePlot=FALSE
waterName='water'
fluc.comp=TRUE
waterMean=1100

test_that('makeSubstrateMat',{
    expect_error(makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag=0,RS.frac,makePlot,saveFig=FALSE,waterName,waterMean))

    expect_equal(ncol(makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag,RS.frac,makePlot,saveFig=FALSE,waterName,waterMean)),length(substrateNames)+2)

#check length is within +/-1 of num.days*24/time.step.h
    expect_lte(length(makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag,RS.frac,makePlot,saveFig=FALSE,waterName,waterMean)[,1]),(num.days*24/time.step.h)+1)

    expect_gte(length(makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag,RS.frac,makePlot,saveFig=FALSE,waterName,waterMean)[,1]),(num.days*24/time.step.h)-1)

#check names
    expect_named(makeSubstrateMat(substrateNames,substrateMeans,num.days,breakfast.start,lunch.start,dinner.start,time.step.h,meal.duration.h,fluc.comp,time.to.reach.colon.h,gamma.mag,RS.frac,makePlot,saveFig=FALSE,waterName,waterMean)[1,],c('time',substrateNames,waterName))
})
