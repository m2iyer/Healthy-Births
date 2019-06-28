## read in dataset

chds_births <- read.csv('chds_births2.csv')


library(plyr)
impute.med <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
no_na <- sapply(chds_births, function(x) { if(is.numeric(x)) impute.med(x) else x }
)

no_na = data.frame(no_na)

na_count = apply(no_na, 2, function(x) any(is.na(x)))

no_na$meth[no_na$meth == 1] <- 0
no_na$meth[no_na$meth == 2] <- 0
no_na$meth[no_na$meth == 3] <- 0
no_na$meth[no_na$meth == 4] <- 0
no_na$meth[no_na$meth == 5] <- 0

no_na$feth[no_na$feth == 1] <- 0
no_na$feth[no_na$feth == 2] <- 0
no_na$feth[no_na$feth == 3] <- 0
no_na$feth[no_na$feth == 4] <- 0
no_na$feth[no_na$feth == 5] <- 0

no_na$feth <- as.factor(no_na$feth)
no_na$meth <- as.factor(no_na$meth)

## Change mht and mwt to BMI
no_na$bmi <- (no_na$mwt/2.2)/((no_na$mht*2.54/100)^2)

## allyssas code

## change categoricals to words
timecats <-c("never","still smokes","during pregnancy","less than year","1-2 yrs","2-3 yrs","3-4 yrs","5-9 yrs",">10 yrs","quit unknown")


no_na$time_cat <- timecats[no_na$time+1]

## first model
chds_births <- no_na[c(1:5, 9, 12, 13, 17:18)]
MQ1 <- lm(wt~(. - marital - income)^2 - feth:bmi - 
            meth:feth + marital + income, data=chds_births)

MModel <- MQ1$model

M1 <- lm(wt ~ 0 + time_cat + gestation + time_cat:gestation, data = chds_births)
summary(M1)

MQ <- lm(wt~(. - marital - income - time_cat)^2 - feth:bmi - 
           meth:feth + marital + income + time_cat + time_cat:gestation, data=chds_births)

M0 <- lm(wt~1, data = chds_births)
#MQInit <- lm(wt~(.)^2, data=chdbirth3)
beta.max <- coef(MQ)
names(beta.max)[is.na(beta.max)]
anyNA(coef(MQ))
MStart <- lm(wt~., data=chds_births)
# time and smoke were taken out
system.time({
  Mfwd <- step(object =M0,
               scope = list(lower = M0, upper = MQ),
               direction = "forward",
               trace = FALSE)
  
})

system.time({
  Mbck <- step(object =MQ,
               scope = list(lower = M0, upper = MQ),
               direction = "backward",
               trace = FALSE)
  
})
system.time({
  Mstp <- step(object =MStart,
               scope = list(lower = M0, upper = MQ),
               direction = "both",
               trace = FALSE)
  
})

## cross validation


# models to compare
M2 <- Mstp
Mnames <- expression(M[TIME], M[STEP])
# Cross-validation setup
nreps <- 2e3 # number of replications
ntot <- nrow(chds_births) # total number of observations
ntrain <- 1000 # size of training set
ntest <- ntot-ntrain # size of test set
mspe1 <- rep(NA, nreps) # sum-of-square errors for each CV replication
mspe2 <- rep(NA, nreps)
logLambda <- rep(NA, nreps) # log-likelihod ratio statistic for each replication
system.time({
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    # refit the models on the subset of training data; ?update for details!
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    # out-of-sample residuals for both models
    # that is, testing data - predictions with training parameter
    M1.res <- chds_births$wt[-train.ind] -
      predict(M1.cv, newdata = chds_births[-train.ind,])
    M2.res <- chds_births$wt[-train.ind] -
      predict(M2.cv, newdata = chds_births[-train.ind,])
    # mean-square prediction errors
    mspe1[ii] <- mean(M1.res^2)
    mspe2[ii] <- mean(M2.res^2)
    # out-of-sample likelihood ratio
    M1.sigma <- sqrt(sum(resid(M1.cv)^2)/ntest) # MLE of sigma
    M2.sigma <- sqrt(sum(resid(M2.cv)^2)/ntest)
    # since res = y - pred, dnorm(y, pred, sd) = dnorm(res, 0, sd)
    logLambda[ii] <- sum(dnorm(M1.res, mean = 0, sd = M1.sigma, log = TRUE))
    logLambda[ii] <- logLambda[ii] -
      sum(dnorm(M2.res, mean = 0, sd = M2.sigma, log = TRUE))
  }
})

# plot rMSPE and out-of-sample log(Lambda)
par(mfrow = c(1,2))
par(mar = c(4.5, 4.5, .1, .1))
boxplot(x = list(sqrt(mspe1), sqrt(mspe2)), names = Mnames, cex = .7,
        ylab = expression(sqrt(MSPE)), col = c("yellow", "orange"))
hist(logLambda, breaks = 50, freq = FALSE,
     xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(logLambda), col = "red") # average value