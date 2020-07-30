#visits_prop.df <- read.csv('./data/visits_Conference.owl_prop_3.csv')
#visits_ant_prop.df <- read.csv('./data/visits_ant_Conference.owl_prop_3.csv')
#visits_order_prop.df <- read.csv('./data/visits_order_Conference.owl_prop_3.csv')

visits_class.df <- read.csv('./data/visits_Conference.owl_class_500.csv')
visits_ant_class.df <- read.csv('./data/visits_ant_Conference.owl_class_500.csv')
visits_order_class.df <- read.csv('./data/visits_order_Conference.owl_class_500.csv')

head(visits_class.df)
head(visits_ant_class.df)
head(visits_order_class.df)

library(coda) 
library(knitr) 
library(kableExtra) 
library(faraway) 
library(mvtnorm)
suppressPackageStartupMessages(
library(rstan,quietly=TRUE, warn.conflicts = FALSE, verbose = FALSE))
#data(prostate, package="faraway") 
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE) 
set.seed(8675309)

getOption("contrasts")
#n <- length(visits.df$visits) 
#k <- 2
#x <- 17 # time to predict
#node_num <- as.numeric(factor(visits_class.df$node))
#visits.df$node_num <- node_num

#aggregate to generate avg hops per node (df1), visits per loop (df3) and avg hops per loop per node (df4)

########################
visits_orde_c_temp.df <- subset(visits_order_class.df, select = -c(nodeType,nodeLevel,nodeName ))
head(visits_orde_c_temp.df)

agg = aggregate(visits_orde_c_temp.df,
                by = list(visits_orde_c_temp.df$nodeId),
                FUN = mean)

typeof(agg)
agg$hop

visits_orde_c_temp1.df <- subset(visits_order_class.df, select = -c(nodeType,nodeLevel,nodeName,hop ))
#visits_orde_c_temp1.df <- subset(visits_order_class.df, select = -c(nodeType,nodeLevel,hop ))
head(visits_orde_c_temp1.df)

aggAvgLoop = aggregate(visits_orde_c_temp.df,
                by = list(visits_orde_c_temp.df$nodeId, visits_orde_c_temp.df$loop),
                FUN = mean)

aggAvgLoop

df4 <- as.data.frame(aggAvgLoop)

colnames(df4)

df4 <- subset(df4, select = -c(Group.1,Group.2) )
names(df4)[1] <- "nodeNumber"


colnames(df4)
head(df4)

nr.of.appearances <- aggregate(x = visits_orde_c_temp1.df, 
                               by = list(unique.values = visits_orde_c_temp1.df$nodeId), 
                               FUN = length)

nr.of.appearances


nr.of.appearancesLoop <- aggregate(x = visits_orde_c_temp1.df, 
                               by = list(visits_orde_c_temp1.df$nodeId, visits_orde_c_temp1.df$loop), 
                               FUN = length)

nr.of.appearancesLoop

df3 <- as.data.frame(nr.of.appearancesLoop)

colnames(df3)

df3 <- subset(df3, select = -c(loop) )
names(df3)[1] <- "nodeNumber"
names(df3)[2] <- "loop"
names(df3)[3] <- "visitsLoop"

colnames(df3)
head(df3)

#visits_class3.df<-merge(x=df3,y=df4,by="nodeNumber","loop",all=FALSE)

df3$hop <- df4$hop 

head(df3)

#visits_class.df$hopAvg <- agg$hop

#do.call(rbind, lapply(agg[,2:4], as.data.frame))
#df <- data.frame(matrix(unlist(ls), ncol = max(lengths(ls)), byrow = TRUE))

# to use with vsits total count
df1 <- as.data.frame(agg)

colnames(df1)

names(df1)[2] <- "nodeNumber"

df1 <- subset(df1, select = -c(loop,Group.1) )

colnames(df1)

#visits_class.df<-merge(x=visits_class.df,y=df1,by="nodeNumber",all.x=TRUE)
#visits_class2.df<-merge(x=visits_class.df,y=df1,by="nodeNumber",all.x=TRUE)




df5 <- merge(x=df3, y=visits_class.df,by="nodeNumber",all.x=TRUE)

visits_class.df <- df5

head(visits_class.df)

names(df3)[3] <- "visits"

head(df3)

#visits_class.df <- df3
###############################

# now starts Gibbs !!! 
head(visits_class.df)
tail(visits_class.df)

par(mfrow=c(1, 1))
plot(visits ~ nodeNumber, visits_class.df,ylab="number vistss")

stripchart(hop ~ nodeNumber, visits_class.df, vertical=TRUE,
           method="stack",xlab="node",ylab="visits")

print(is.factor(visits_class.df$nodeNumber))

visits_class.df$nodeNumber <- factor(visits_class.df$nodeNumber)

print(is.factor(visits_class.df$nodeNumber))

# it is not working with you have NAs in the summary
visits_c.lm <- lm(visits~node + loop + hop, data=visits_class.df)
sumvisits_c.lm <- summary(visits_c.lm)
sumvisits_c.lm

#cell mean model only to check values 
visits_c_cellmeanmodel.lm <- lm(visits~nodeNumber + loop + hop -1  , data=visits_class.df)
summary(visits_c_cellmeanmodel.lm)

#visits_p.lm <- lm(visits~node, data=visits_prop.df)
#summary(visits_p.lm)

#class
X <- model.matrix(visits_c.lm) 
y <- visits_c.lm$model[,1]
n <- nrow(X)
p <- ncol(X)

XtX<- t(X) %*% X
bhat<- coef(visits_c.lm) 
ybar<- mean(visits_class.df$visits) 
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

parms <- rownames(sumvisits_c.lm$coefficients)

par(mfrow=c(2, 3),mar = rep(2, 4))

# for until nrow(W)
nrow(W)
for(j in 1:nrow(W)){
  #plot(1:M, beta1[j,-1], type="l", xlab="", ylab="", main = colnames(visits_class.df)[j])  
  plot(1:M, beta1[j,-1], type="l", xlab="", ylab="", main = colnames(parms[j]))
  lines(1:M, beta2[j,-1], col = 2, lty = 2)
  lines(1:M, beta3[j,-1], col = 4, lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}; for(j in 1:nrow(W)){
  plot(density(na.omit(beta1[j,])), main = NA)
  lines(density(na.omit(beta2[j,])), col = 2, lty = 2)
  lines(density(na.omit(beta3[j,])), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}

stats <- c('mean', 'med', 'sd', 'x2.5%', 'x97.5%')
niceTables <- rbind(
  data.frame(grp = 'lm', stats = c('mean', 'sd'), 
             matrix(t(summary(visits_c.lm)$coefficients[, c('Estimate', 'Std. Error')]), nrow = 2)),
  data.frame(grp = 'Gibbs', stats, apply(
    cbind(beta1[, -1], beta2[, -1], beta3[, -1]),
    1, function(x) {
      c(mean(x), median(x), sd(x), quantile(x, c(.025, .975)))})))
#colnames(niceTables) <- c('gr', 'stats', colnames(visits_class.df))
colnames(niceTables) <- c('gr', 'stats','nodeNumber','loop','visits','hop')

kable(xtabs(visits ~ gr + stats, data = niceTables))
kable(xtabs(nodeNumber ~ gr + stats, data = niceTables))
kable(xtabs(hop ~ gr + stats, data = niceTables))
kable(xtabs(loop ~ gr + stats, data = niceTables))

######
# MCMC 

zvisits.stanc <- stanc(file="./Stan/visits.stan")
zvisits.stanmod <- stan_model(stanc_ret=zvisits.stanc)

zvisits.data <- list(
  N=dim(X)[1],
  p=dim(X)[2],
  X=X,
  y=visits_c.lm$model[,1],
  m0=rep(0, dim(X)[2]),
  prec0=diag(rep(1/V, ncol(X))),
  nu0=1/1000,
  sig20 = 1)

set.seed(5551212)
binits1 <- rmvnorm(n=1,mean=coef(visits_c.lm), sigma=vcov(visits_c.lm)/4)
binits2 <- rmvnorm(n=1,mean=coef(visits_c.lm), sigma=vcov(visits_c.lm))
binits3 <- rmvnorm(n=1,mean=coef(visits_c.lm), sigma=vcov(visits_c.lm)*4)
attach(zvisits.data)
sigma2inits1 <- sum((y - X%*%binits1[1,])^2)/(N-p)
sigma2inits2 <- sum((y - X%*%binits2[1,])^2)/(N-p)
sigma2inits3 <- sum((y - X%*%binits3[1,])^2)/(N-p)
detach(zvisits.data)

zvisits.init <- list(
  list(beta=binits1[1,], sigma2 = sigma2inits1),
  list(beta=binits2[1,], sigma2 = sigma2inits2),
  list(beta=binits3[1,], sigma2 = sigma2inits3)
)

zvisits.fit <- rstan::sampling(zvisits.stanmod,
                                 seed=8675309, warmup = 5000,
                                 data=zvisits.data,
                                 pars=c("beta","sigma2"),
                                 chains = 3, iter=10000,
                                 init = zvisits.init,
                                 refresh = 1000)

(zvisits.sum<- summary(zvisits.fit))
mcmc <- As.mcmc.list(zvisits.fit)

nchain(mcmc)
niter(mcmc)
nvar(mcmc)
varnames(mcmc)

for(j in 1:nvar(mcmc)){
  plot(as.numeric(mcmc[[1]][, j]), type="l", xlab="", ylab="", main = names(zvisits.fit)[j])
  lines( mcmc[[2]][, j], col = alpha(2, .5), lty = 2)
  lines( mcmc[[3]][, j], col = alpha(4, .5), lty = 2)
  legend('topleft', legend = expression(paste(beta[j], " | y")), bty = 'n')
}; for(j in 1:nvar(mcmc)){
  plot(density(mcmc[[1]][, j]), main = NA)
  lines(density(mcmc[[2]][, j]), col = 2, lty = 2)
  lines(density(mcmc[[3]][, j]), col = 4, lty = 4)
  legend('topright',legend = c(expression(sigma/4 ), expression(sigma), expression(sigma * 4)),
         bty = 'n', lty = 1:3, col = c(1:2, 4))
}

# trying another plot

visits.array<- rstan::extract(zvisits.fit, permuted=FALSE, inc_warmup=TRUE)
visits.mcmc.list<- vector("list", 3)
for(chain in 1:3){
  visits.mcmc.list[[chain]]<- coda::as.mcmc(visits.array[,chain,])
}
visits.mcmc.list <- coda::as.mcmc.list(visits.mcmc.list)
coda::niter(visits.mcmc.list)
par(mar=c(2,2,2,2))
coda:::plot.mcmc.list(visits.mcmc.list)
coda:::summary.mcmc.list(window(visits.mcmc.list, start=5001, end=10000))
