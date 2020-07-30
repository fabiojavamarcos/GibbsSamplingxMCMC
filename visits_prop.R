visits_prop.df <- read.csv('./data/visits_Conference.owl_prop_3.csv')
visits_ant_prop.df <- read.csv('./data/visits_ant_Conference.owl_prop_3.csv')
visits_order_prop.df <- read.csv('./data/visits_order_Conference.owl_prop_3.csv')

visits_class.df <- read.csv('./data/visits_Conference.owl_class_3.csv')
visits_ant_class.df <- read.csv('./data/visits_ant_Conference.owl_class_3.csv')
visits_order_class.df <- read.csv('./data/visits_order_Conference.owl_class_3.csv')

visits_prop.df
visits_ant_prop.df
visits_order_prop.df

library(coda) 
library(knitr) 
library(kableExtra) 
library(faraway) 
library(mvtnorm)
suppressPackageStartupMessages(
  library(rstan,quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE))
#data(prostate, package="faraway") 
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE) set.seed(8675309)


#n <- length(visits.df$visits) 
#k <- 2
#x <- 17 # time to predict
#node_num <- as.numeric(factor(visits_class.df$node))
#visits.df$node_num <- node_num

visits_c.lm <- lm(visits~nodeNumber, data=visits_class.df)
summary(visits_c.lm)

visits_p.lm <- lm(visits~nodeNumber, data=visits_prop.df)
summary(visits_p.lm)

#class
X <- model.matrix(visits_p.lm) 
y <- visits_c.lm$model[,1]
n <- nrow(X)
p <- ncol(X)

XtX<- t(X) %*% X
bhat<- coef(visits.lm) 
ybar<- mean(visits.df$visits) 
nu0 <- ybar
yhat<- X %*%bhat 
epshat<-y-yhat 
sigmahat<-sum(epshat^2)/(n-p)

sig02<- 0; ## prior scale squared 
m0<- rep(0,p) ## prior beta mean 
V<- (log(10)/1.96)^2
Sigma0inv<- diag(c(0,rep(1/V,p-1))) 
nuhat<- nu0 + n
M <- 10000 ## MCMC iterations
beta_all <- list()
sigma2_all <- matrix(rep(NA, (M+1) * 3), ncol = 3) 
sigma2_all[1, 1]<- sigmahat^2/4 ## initial sigma20 
sigma2_all[1, 2]<- sigmahat^2 ## initial sigma20 
sigma2_all[1, 3]<- 4 * sigmahat^2 ## initial sigma20 
sXtX <- solve(XtX)

#visits.list<- list(n=n, k=k, x=x, node=visits.df$node, visits=visits.df$visits)

#writeLines(readLines("visits.stan"))


#nchains <- 3 
#visits.data.list <- list(
# n = visits.list$n,
#k = visits.list$k,
#x = visits.list$x,
#node = visits.list$node, 
#visits = visits.list$visits
#)

for (s in 1:3){
  sigma2 <- as.numeric(sigma2_all[, s]) 
  beta <- matrix(NA,nrow=p,ncol=M+1) 
  for (i in 1:M){
    W <- solve(sigma2[i] * Sigma0inv + XtX)%*%XtX
    mhat <- (diag(p) - W) %*% m0 + W %*% bhat #
    Sighat <- W %*% sXtX
    beta[, i+1] <- as.vector(rmvnorm(n=1,mean=mhat, sigma=sigma2[i]*Sighat)) 
    sigma2hat <- (nu0*sig02 + sum((y - X%*%beta[,i+1])^2))/nuhat
    sigma2[i+1] <- nuhat * sigma2hat / rchisq(n=1, df=nuhat) }
  beta_all[[s]] <- beta
  beta_all[[s]] <- beta
  sigma2_all[, s] <- sigma2
}
sigma2.1 <- sigma2_all[, 1]
sigma2.2 <- sigma2_all[, 2]
sigma2.3 <- sigma2_all[, 3]
beta1 <- beta_all[[1]]
beta2 <- beta_all[[2]]
beta3 <- beta_all[[3]]

par(mfrow=c(2, 3),mar = rep(2, 4))
for(j in 1:2){
  plot(1:M, beta1[j,-1], type="l", xlab="", ylab="", main = colnames(visits.df)[j])
  lines(1:M, beta2[j,-1], col = 2, lty = 2)
  lines(1:M, beta3[j,-1], col = 4, lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}; for(j in 1:2){
  plot(density(na.omit(beta1[j,])), main = NA)
  lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
  lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}




