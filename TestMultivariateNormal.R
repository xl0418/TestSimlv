#library(devtools)
#install_github("xl0418/SimLV")


library(SimLV)
library(ggplot2)
library(reshape2)


# Simulate data
# test
n0 <- c(1, 1, 1)
r <- c(1, 1, 1)
alpha <- matrix(c(0.2, 0.1, 0.3, 0.1, 0.5, 0.3, 0.1, 0.12, 0.3), nrow = 3, byrow = TRUE)
tmax <- 100
simresult <- simlv(n0 = n0, r = r, alphaij = alpha, tmax = tmax)

# long data frame
df <- melt(simresult, id.vars=c("time"))

# plot the result
ggplot(df, aes(x = time, y = value, group = variable, color = factor(variable))) + geom_line(linewidth = 2) + theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species")




library(mgcv)
library(dplyr)

# add some lag values to the data
simresult<-simresult%>%mutate(X1lag=lag(X1),X2lag=lag(X2),X3lag=lag(X3),X1diff=X1-lag(X1),X2diff=X2-lag(X2),X3diff=X3-lag(X3))
simresult<-simresult[complete.cases(simresult),]

#first fit intercept only , no lags allowed
f <- list(X1~1,X2~1,X3~1)
b1 <- gam(f, family = mvn(d = 3), data = simresult)

predictb1<-predict(b1,type="response") # prediction
b1_covariance<-solve(crossprod(b1$family$data$R)) #covariance matrix

##compare alpha and covariance matrix
L1<-data.frame(processmodel=c(alpha),mvn_fit=c(b1_covariance),type=c("1,1","1,2","1,3","1,2","2,2","2,3","1,3","2,3","3,3"))

ggplot(L1, aes(x=mvn_fit,y=processmodel,color=type))+geom_point()+geom_smooth(color="black",method="lm")

#gives the overall marginal means
ggplot(data.frame(time=simresult$time,values=c(simresult$X1,simresult$X2,simresult$X3,predictb1[,1],predictb1[,2],predictb1[,3]),species=rep(c("X1","X2","X3"),each=nrow(simresult)),type=rep(c("observed","fitted"),each=3*nrow(simresult))), aes(x = time, y = values,linetype=factor(type), color = factor(species))) + geom_line(linewidth = 2) + theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species")

# now a model with AR1 autocorrelation
f1 <- list(X1~X1lag+1,X2~X2lag+1,X3~X3lag+1)
b2 <- gam(f1, family = mvn(d = 3), data = simresult)

predictb2<-predict(b2,type="response")
b2_covariance<-solve(crossprod(b2$family$data$R))

L2<-data.frame(processmodel=c(alpha),mvn_fit=c(b2_covariance),type=c("1,1","1,2","1,3","1,2","2,2","2,3","1,3","2,3","3,3"))

ggplot(L2, aes(x=mvn_fit,y=processmodel,color=type))+geom_point()+geom_smooth(color="black",method="lm")

ggplot(data.frame(time=simresult$time,values=c(simresult$X1,simresult$X2,simresult$X3,predictb2[,1],predictb2[,2],predictb2[,3]),species=rep(c("X1","X2","X3"),each=nrow(simresult)),type=rep(c("observed","fitted"),each=3*nrow(simresult))), aes(x = time, y = values,linetype=factor(type), color = factor(species))) + geom_line(linewidth = 2) + theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species")

# now a model including first lags of all species
f2 <- list(X1~X1lag+X2lag+X3lag+1,X2~X1lag+X2lag+X3lag+1,X3~X1lag+X2lag+X3lag+1)
b3 <- gam(f2, family = mvn(d = 3), data = simresult)

predictb3<-predict(b3,type="response")
b3_covariance<-solve(crossprod(b3$family$data$R))

L3<-data.frame(processmodel=c(alpha),mvn_fit=c(b3_covariance),type=c("1,1","1,2","1,3","1,2","2,2","2,3","1,3","2,3","3,3"))

ggplot(L3, aes(x=mvn_fit,y=processmodel,color=type))+geom_point()+geom_smooth(color="black",method="lm")

ggplot(data.frame(time=simresult$time,values=c(simresult$X1,simresult$X2,simresult$X3,predictb3[,1],predictb3[,2],predictb3[,3]),species=rep(c("X1","X2","X3"),each=nrow(simresult)),type=rep(c("observed","fitted"),each=3*nrow(simresult))), aes(x = time, y = values,linetype=factor(type), color = factor(species))) + geom_line(linewidth = 2) + theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species")

coefs<-coef(b3)[c("X1lag","X2lag","X3lag","X1lag.1","X2lag.1","X3lag.1","X1lag.2","X2lag.2","X3lag.2")]
#coefs<-coef(b3)[c("X1lag","X1lag.1","X1lag.2","X2lag","X2lag.1","X2lag.2","X3lag","X3lag.1","X3lag.2")]

L3x<-data.frame(processmodel=c(alpha),mvn_fit=coefs,type=c("1,1","1,2","1,3","1,2","2,2","2,3","1,3","2,3","3,3"))
ggplot(L3x, aes(x=mvn_fit,y=processmodel,color=type))+geom_point()+geom_smooth(color="black",method="lm")

# I guess this gives good fits, but is not the way to go to get the alphas correct.
#when including lagged abundances the covariances go to zero, but maybe not when the error is multiplicative in the generation (see comment below)
summary(b1)
summary(b2)
summary(b3)


## Maybe we need to probably model first differences of abundances instead of abundances
# and find a way to make the process model additive instead of multiplicative? Can log transforms help here?
# also the error in your function is currently additive. 
# this 1) can lead to negative values and 
# 2) probably in nature it is rather multiplicative? say proportional to the true value rather than absolute?
# would a multiplicative error also make things easier when log transform things (not sure how we tackle negative dn/dt) 
# or we could just model the multiplicator per time step? say its N+N*r*(1-bla) but could formulate as N*(1+N^2*r*(1-bla)) where (1+.. is) always positive?
# then do -1 and take log so log(N^2*r*(1-bla)) which is 2*log(N)+log(r)+log(1-bla)
# maybe nothing makes sense and thats totally nonsense :D just writing while thinking , lets chat on slack. 

