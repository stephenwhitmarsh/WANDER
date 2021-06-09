install.packages("lme4")
install.packages("mediation")
install.packages("mediation", dependencies=TRUE)
install.packages("mvtnorm")
install.packages("lmerTest")
install.packages("ggplot2")
install.packages("MuMIn") # for REVISION

# data <- read.csv(file="C:/Users/Admin/Dropbox/WANDER_article/trialdata_triallength.csv", sep=';', header=TRUE, dec=',')
# data <- read.csv(file="C:/Users/Admin/Dropbox/WANDER_article/trialdata_triallength_2.csv", sep=';', header=TRUE, dec=',')
# data <- read.csv(file="Z:/WANDER/data/trialdata_triallength_2.csv", sep=',', header=TRUE, dec='.')
# data <- read.csv(file="D:/Dropbox/WANDER_article/trialdata_triallength_3.csv", sep=',', header=TRUE, dec='.')
data <- read.csv(file="C:/Users/Stephen Whitmarsh/Dropbox/WANDER_article/trialdata_triallength_3.csv", sep=',', header=TRUE, dec='.')
data <- read.csv(file="C:/Users/Stephen Whitmarsh/Dropbox/WANDER_article/trialdata_triallength_3.csv", sep=',', header=TRUE, dec='.')
data <- read.csv(file="D:/Dropbox/ENS/WANDER_article/trialdata_triallength_3.csv", sep=',', header=TRUE, dec='.')



library(lme4)
library(ggplot2)

data$suj    <- factor(data$SubjectID)
data$sS     <- data$SteadyState # because jerome doesnt like long names
data$Rating <- 9 - data$Rating
par(mfrow=c(2,2))
sujs <- unique(data$suj)

# correct time as minutes to get into more of the same range as rest of predictors
data$corTime <- data$Timestamp / 600000

# We z-score by hand
e <- evalq(aggregate(list(meanP=Pupil), list(suj=suj), mean), data)
f <- evalq(aggregate(list(sdP=Pupil), list(suj=suj), sd), data)
data <- merge(data, e)
data <- merge(data, f)

data$zPupil <- (data$Pupil-data$meanP)/data$sdP
e <- evalq(aggregate(list(meanR=Rating), list(suj=suj), mean), data)
f <- evalq(aggregate(list(sdR=Rating), list(suj=suj), sd), data)

data <- merge(data, e)
data <- merge(data, f)
data$zRating <- (data$Rating-data$meanR)/data$sdR

e <- evalq(aggregate(list(meanA=Alpha), list(suj=suj), mean), data)
f <- evalq(aggregate(list(sdA=Alpha), list(suj=suj), sd), data)
data <- merge(data, e)
data <- merge(data, f)

data$zAlpha <- (data$Alpha-data$meanA)/data$sdA
e <- evalq(aggregate(list(meanS=sS), list(suj=suj), mean), data)
f <- evalq(aggregate(list(sdS=sS), list(suj=suj), sd), data)

data <- merge(data, e)
data <- merge(data, f)
data$zSteady <- (data$sS-data$meanS)/data$sdS

e <- evalq(aggregate(list(meanRMS=RMS), list(suj=suj), mean), data)
f <- evalq(aggregate(list(sdRMS=RMS), list(suj=suj), sd), data)
data <- merge(data, e)
data <- merge(data, f)
data$zRMS <- (data$RMS-data$meanS)/data$sdS

# add ratings, pupil and alpha of previous trial
data$pZRating  <- c(NA, data$zRating)[1:nrow(data)]
data$pZPupil   <- c(NA, data$zPupil)[1:nrow(data)]
data$pZAlpha   <- c(NA, data$zAlpha)[1:nrow(data)]
data$pZSteady   <- c(NA, data$zSteady)[1:nrow(data)]

data$pZRMS     <- c(NA, data$zRMS)[1:nrow(data)]
data$pZRating2 <- c(NA, data$pZRating)[1:nrow(data)]
data$pZPupil2  <- c(NA, data$pZPupil)[1:nrow(data)]
data$pZAlpha2  <- c(NA, data$pZAlpha)[1:nrow(data)]
data$pZRMS2    <- c(NA, data$pZRMS)[1:nrow(data)]
data$pZRating3 <- c(NA, data$pZRating2)[1:nrow(data)]
data$pZPupil3  <- c(NA, data$pZPupil2)[1:nrow(data)]
data$pZAlpha3  <- c(NA, data$pZAlpha2)[1:nrow(data)]
data$pZRMS3    <- c(NA, data$pZRMS2)[1:nrow(data)]
data$pZRating4 <- c(NA, data$pZRating3)[1:nrow(data)]
data$pZPupil4  <- c(NA, data$pZPupil3)[1:nrow(data)]
data$pZAlpha4  <- c(NA, data$pZAlpha3)[1:nrow(data)]
data$pZRMS4    <- c(NA, data$pZRMS3)[1:nrow(data)]
data <- droplevels(data[!is.na(data$pZAlpha4),])
head(data)


## Correlation of predictors

library(lmerTest)

# de the tests of ratings vs. time and block previously done in MATLAB
l <- lmer(Rating ~ Trialnr + Blocknr + (1|suj), data); summary(l)  # 
l <- lm(Rating ~ Trialnr + Blocknr, data); summary(l)  # 
ggplot(data, aes(x=Trialnr, y=Rating)) + geom_bar(stat = "summary", fun.y = "mean") + facet_wrap(~Blocknr, scales = "free_y", nrow = 1, ncol = 4)
evalq(pairs(~zSteady+zPupil+zAlpha+zRating), data)

# REVISION - USED TO BE DONE IN MATLAB
library(MuMIn)
l1 <- lmer(zRating ~ Length + Trialnr + Blocknr + (1 | suj), data);  summary(l)  
r.squaredGLMM(l1)


# plot correlations in multipanel
# 
# install.packages("rlang")
# install.packages("devtools")
# library(devtools)
# install_github("easyGgplot2", "kassambara")
# library(easyGgplot2)
# 
# plot1 <- ggplot(data, aes(x=zSteady, y=zAlpha)) + geom_point() + geom_smooth(method = "lm")
# plot2 <- ggplot(data, aes(x=zRating, y=zAlpha)) + geom_point() + geom_smooth(method = "lm")
# plot3 <- ggplot(data, aes(x=zRating, y=zSteady)) + geom_point() + geom_smooth(method = "lm")
# 
# ggplot2.multiplot(plot1,plot2,plot3, cols=2)



require(plyr)
func <- function(xx)
{
  return(data.frame(COR = cor(xx$zAlpha, xx$zSteady)))
}

c <- ddply(data, .(suj), func)

ggplot(data = c, aes(x = "", y = COR)) +   geom_boxplot() +
geom_dotplot(binaxis='y', stackdir='center') + ylab('Correaltion Alpha x SteadyState') + xlab('')


# explore relationship between variables
l <- lmer(zSteady ~ zPupil  +(1|suj), data); summary(l)   # pupil does not predict SS
l <- lmer(zAlpha  ~ zPupil  +(1|suj), data); summary(l)   # pupil predict alpha strongly
l <- lmer(zPupil  ~ zAlpha  +(1|suj), data); summary(l)   # alpha predicts pupil strongly
l <- lmer(zSteady ~ zAlpha  +(1|suj), data); summary(l)   # alpha predicts SS very strongly
l <- lmer(zPupil  ~ zSteady  +(1|suj), data); summary(l)  # SS does not predicts pupil
l <- lmer(zAlpha  ~ zSteady  +(1|suj), data); summary(l)  # SS predicts alpha very strongly

# while alpha predicts SS, and they predict rating, SS does not predict rating
l1 <- lmer(zRating ~ zPupil + zAlpha + zSteady + (1 | suj), data); summary(l1) 
l2 <- lmer(zRating ~ zPupil + zAlpha + (1 | suj), data);  summary(l2)  

# model fits better without steady-state
anova(l1,l2)

# explore single predictors
l <- lmer(zRating ~ zAlpha + (1 | suj), data);  summary(l) # yes
l <- lmer(zRating ~ zPupil + (1 | suj), data);  summary(l) # yes
l <- lmer(zRating ~ zSteady + (1 | suj), data); summary(l) # no

# get some p-values
detach(package:lmerTest)
library(lmerTest)
l1 <- lmer(zRating ~ zPupil + zAlpha + zSteady + (1 | suj), data); summary(l1)
l2 <- lmer(zRating ~ zPupil + zAlpha + (1 | suj), data); summary(l2)
anova(l1, l2)

## Sanity check: test model in two steps: first per subject, then t-test. This shows similar results 
L <- lapply(sujs, function (s){
                l <- lm(zRating ~ zPupil + zAlpha + zSteady, data=data[data$suj==s,])
                return(data.frame(suj=s, zPupil=as.numeric(l$coeff[2]), zAlpha=as.numeric(l$coeff[3]),
                                  zSteady=as.numeric(l$coeff[4])))
            })

coefDF <- do.call(rbind, L)
t.test(coefDF$zPupil)
t.test(coefDF$zAlpha)

## test mixed model, now including blocknr and time, with and without steady-state. Similarly better without:  8272.8 - 8280.6 = -7.8; exp(-.5*-7.8) = 49.40245
l0 <- glm(zRating ~ 1, data = data, family = gaussian); summary(l0)
l1 <- lmer(zRating ~ zPupil + zAlpha + zSteady + Blocknr + corTime + (1 | suj), data); summary(l1)
l2 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + (1 | suj), data); summary(l2)
anova(l1, l2)

# REVISION
r.squaredGLMM(l1)
r.squaredGLMM(l2)



## so do it without steady-state, and compare with and without interaction between pupil and alpha (there is none):  8466.3 - 8473.7 = -7.4; exp(-.5*-7.4) = 40.4473
l1 <- lmer(zRating ~ zPupil * zAlpha + Blocknr + corTime + (1 | suj), data); summary(l1)
l2 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + (1 | suj), data); summary(l2)
anova(l1, l2)

## test for interaction block * time; exp(-.5*-7.3) = 38.47, going for non-interaction???? (CHECK WITH JEROME, seems otherwise)
## 8466.3 - 8474.1 = -7.8; exp(-.5*-7.8) = 49.40245
l1 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + (1  | suj), data); summary(l1) # BIC 8272.8; later: 8466.3
l2 <- lmer(zRating ~ zPupil + zAlpha + Blocknr * corTime + (1  | suj), data); summary(l2) # BIC 8280.1l later: 8474.1
detach(package:lmerTest)
anova(l1, l2)

## check mediation http://dspace.mit.edu/handle/1721.1/91154
## Simple explanation and comparison with different methods: http://www.psych.mcgill.ca/perpg/fac/falk/mediation.html#CIcalculator
# data("framing", package = "mediation")
# med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
# out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income, data = framing, family = binomial("probit"))
# med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", robustSE = TRUE, sims = 100)
# # summary(med.out)
# 

library(mediation)
# med.fit <- lmer(zPupil  ~ zAlpha + Blocknr + corTime + (1  | suj), data); summary(med.fit);
# out.fit <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + (1  | suj), data); summary(out.fit);
# med.out <- mediate(med.fit, out.fit, treat = "zRating", mediator = "zPupil")
# https://stackoverflow.com/questions/12485309/r-mediation-analysis-bootstrapping
##

# needed for error: https://stackoverflow.com/questions/41253257/error-using-mediation-package-with-lme4-model-mediator-model-is-not-yet-impleme
detach(package:lmerTest)
library(lme4) 

# model ratings as response: effect of pupil on rating, modulated by alpha
med.fit <- lmer(zAlpha  ~ zPupil + Blocknr + corTime + (1  | suj), data = data); summary(med.fit);
out.fit <- lmer(zRating ~ zAlpha + zPupil  + Blocknr + corTime + (1  | data$suj), data = data); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);
plot(med.out)

## To calculate With calculator: http://www.psych.mcgill.ca/perpg/fac/falk/mediation.html#CIcalculator
A.fit <- lmer(zAlpha  ~ zPupil + Blocknr + corTime + (1  | suj), data = data); summary(A.fit);
B.fit <- lmer(zRating ~ zAlpha  + zPupil + Blocknr + corTime + (1  | data$suj), data = data); summary(B.fit);
C.fit <- lmer(zRating ~ zPupil  + Blocknr + corTime + (1  | data$suj), data = data); summary(C.fit);

# to calculate DF
library(lmerTest)



# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.00775      0.00311         0.01   0.002 ** 
# ADE             0.19786      0.16410         0.23  <2e-16 ***
# Total Effect    0.20561      0.17150         0.24  <2e-16 ***
# Prop. Mediated  0.03691      0.01522         0.06   0.002 ** 

med.fit <- lmer(zAlpha ~ zPupil + (1  | suj), data = data); summary(med.fit);
out.fit <- lmer(zRating ~ zPupil + zAlpha + (1  | data$suj), data = data); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);

# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            -0.0183      -0.0265        -0.01  <2e-16 ***
# ADE             -0.1290      -0.1632        -0.09  <2e-16 ***
# Total Effect    -0.1473      -0.1805        -0.11  <2e-16 ***
# Prop. Mediated   0.1235       0.0738         0.19  <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



## model steady-state power as response

med.fit <- lmer(zAlpha ~ zPupil +  Blocknr + corTime + (1  | suj), data = data); summary(med.fit);
out.fit <- lmer(zSteady ~ zPupil + zAlpha + Blocknr + corTime + (1  | data$suj), data = data); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);

## To calculate With calculator: http://www.psych.mcgill.ca/perpg/fac/falk/mediation.html#CIcalculator
A.fit <- lmer(zAlpha  ~ zPupil + Blocknr + corTime + (1  | suj), data = data); summary(A.fit);
B.fit <- lmer(zSteady ~ zAlpha  + zPupil + Blocknr + corTime + (1  | data$suj), data = data); summary(B.fit);
C.fit <- lmer(zSteady ~ zPupil  + Blocknr + corTime + (1  | data$suj), data = data); summary(C.fit);

# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.000169    -0.002958         0.00    0.92    
# ADE             0.144238     0.107943         0.18  <2e-16 ***
# Total Effect    0.144407     0.108415         0.18  <2e-16 ***
# Prop. Mediated  0.000914    -0.020231         0.02    0.92    
#---
  
med.fit <- lmer(zPupil ~ zAlpha + (1  | suj), data = data); summary(med.fit);
out.fit <- lmer(zSteady ~ zPupil + zAlpha + (1  | data$suj), data = data); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);

# Estimate 95% CI Lower 95% CI Upper p-value
# ACME            0.00000      0.00000         0.00    1.00
# ADE            -0.00253     -0.03822         0.03    0.86
# Total Effect   -0.00253     -0.03822         0.03    0.86
# Prop. Mediated  0.00000      0.00000         0.00    1.00



# model RMS power as response

data_nonan <- data[!is.na(data$zRMS), ]

med.fit <- lmer(zPupil ~ zAlpha + (1  | suj), data = data_nonan); summary(med.fit);
out.fit <- lmer(zRMS ~ zPupil + zAlpha + (1  | data_nonan$suj), data = data_nonan); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zAlpha", mediator = "zPupil"); summary(med.out);

# Estimate 95% CI Lower 95% CI Upper p-value   
# ACME            -0.0175      -0.0311        -0.01   0.002 **
# ADE             -0.1349      -0.2479        -0.02   0.016 * 
# Total Effect    -0.1524      -0.2602        -0.04   0.008 **
# Prop. Mediated   0.1130       0.0411         0.45   0.010 **

# model steady state power as response
data_nonan <- data[!is.na(data$zSteady), ]

med.fit <- lmer(zPupil ~ zAlpha + (1  | data_nonan$suj), data = data_nonan); summary(med.fit);
out.fit <- lmer(zSteady ~ zPupil + zAlpha + (1  | data_nonan$suj), data = data_nonan); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);
plot(med.out)


# Estimate 95% CI Lower 95% CI Upper p-value
# ACME            0.00000      0.00000         0.00    1.00
# ADE            -0.00251     -0.03799         0.03    0.88
# Total Effect   -0.00251     -0.03799         0.03    0.88
# Prop. Mediated  0.00000      0.00000         0.00    1.00

# med.fit <- lm(data$zAlpha ~ data$zPupil ); summary(med.fit);
# out.fit <- lm(data$zRating ~ data$zAlpha + data$zPupil ); summary(out.fit);
# data_s1 <- data[data$suj == 1,]
# med.fit <- lm(zAlpha ~ zPupil, data = data_s1 ); summary(med.fit);
# out.fit <- lm(zRating ~ zAlpha + zPupil, data = data_s1 ); summary(out.fit);

med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha"); summary(med.out);


# https://stackoverflow.com/questions/48708168/error-in-r-data-framem-data-treat-undefined-columns-selected-runni?rq=1
out.fit <- lm(zRating ~ zPupil, data); summary(out.fit);
med.fit <- lm(zRating ~ zPupil + zAlpha, data); summary(med.fit);
anova(med.fit,out.fit)
BIC(med.fit,out.fit)
AIC(med.fit,out.fit)

med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha", boot=TRUE, sims=500); summary(med.out)

# meditator mixed model not implemented
out.fit <- lmer(zRating ~ zPupil + zAlpha + (1  | suj), data); summary(out.fit);
med.fit <- lmer(zRating ~ zPupil + zAlpha + Blocknr + (1  | suj), data); summary(med.fit);
med.out <- mediate(med.fit, out.fit, treat = "zPupil", mediator = "zAlpha", boot=FALSE, sims=500); summary(med.out)




########

# needed for error: https://stackoverflow.com/questions/41253257/error-using-mediation-package-with-lme4-model-mediator-model-is-not-yet-impleme
detach(package:lmerTest)
library(lme4) 

# model SS as response: effect of alpha on SS, modulated by ratings
med.fit <- lmer(zRating ~ zAlpha + Blocknr + corTime + (1  | suj), data = data); summary(med.fit);
out.fit <- lmer(zSteady ~ zRating + zAlpha  + Blocknr + corTime + (1  | data$suj), data = data); summary(out.fit);
med.out <- mediate(med.fit, out.fit, treat = "zAlpha", mediator = "zRating"); summary(med.out);
plot(med.out)





## TIME LAG ANALYSIS
library(lmerTest)

## add Previous trial 
# 8474.1 - 8403.4 = 70.7; exp(-.5*70.7) = 4.443141e-16
# 8474.1 - 8411.3  = 62.8; exp(-.5*62.8) = 2.307561e-14
l0 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + (1  | suj), data); summary(l0);
l1 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZRating + pZPupil + pZAlpha + (1  | suj), data); summary(l1);
anova(l0, l1)

# REVISION
r.squaredGLMM(l1)


l0 <- lmer(zSteady ~ zPupil + zAlpha + Blocknr + corTime + (1  | suj), data); summary(l0);
l1 <- lmer(zSteady ~ zPupil + zAlpha + Blocknr + corTime + pZSteady + pZPupil + pZAlpha + (1  | suj), data); summary(l1);
anova(l0, l1)

# REVISION
r.squaredGLMM(l1)



l1 <- lmer(zRating ~ Blocknr + corTime + pZRating + pZPupil + pZAlpha + (1  | suj), data); summary(l1);

lpRating <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZRating + (1  | suj), data); summary(lpRating);
lpAlpha1 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZAlpha + (1  | suj), data); summary(lpAlpha1);
lpPupil1 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZPupil + (1  | suj), data); summary(lpPupil1);

lpAlpha2 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZAlpha + pZRating + (1  | suj), data); summary(lpAlpha2);
lpPupil2 <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZPupil + pZRating + (1  | suj), data); summary(lpPupil2); # in this model not only previous ratings, but also previous pupil works out

## 8474.1 - 8411.3 = 62.8; exp(-.5*62.8) = 2.307561e-14
lpAll    <- lmer(zRating ~ zPupil + zAlpha + Blocknr + corTime + pZPupil + pZAlpha  + pZRating + (1  | suj), data); summary(lpAll); # in this model not only previous ratings, but also previous pupil works out
anova(l0, lpAll)
summary(lpAll)


