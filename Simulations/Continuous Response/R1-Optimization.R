
setwd("C:\\Research\\An Alternative to Greedy Search\\Programs\\Simulations\\")
rm(list=ls(all=TRUE))



# ========================================
# FUNCTION rdat GENERATES DATA FOR SPLITS
# ========================================

rdat <- function(beta1=2, cutoff=.5, sd=1, n=100, discrete=F)
{
	if (discrete){
		K <- 50
		x <- sample(1:K,n, replace=T)/K
	} else {
		x <- runif(n)
	}
	y <- beta1*sign(x<=cutoff) + rnorm(n, mean=0, sd=sd)
	y <- y-mean(y)
   	return(data.frame(y=y, x=x))
}



# =====================================
# Plot of the sigmoid function
# =====================================


expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
expit0 <- function(x, c, a=10)  (tanh(a*(x-c)/2)+1)/2 # VERIFIED
# expit <- function(x) exp(x)/(1+exp(x))
x <- -100:100/100
a <- c(5, 10, 20, 30, 50)

postscript(file="fig1.eps", horizontal=F)
par(mfrow=c(1,1), mar=rep(4,4))
plot(c(min(x), max(x)), c(0, 1), type="n", xlab=expression(x), ylab=expression(expit(ax)))
for (i in 1:length(a)){
	a.i <- a[i]
	y <- expit(x*a.i)
	lines(x, y, lty=1, col=i+1)
}
legend(x=-.99, y=.99, legend=paste("a=", a, sep=""), 
	lty=rep(1, 5), lwd=1,
	col=2:(length(a)+1), box.col ="gray84", bg ="gray84")
lines(c(-1, 0, 0, 1), c(-.004, -0.004, 1.004, 1.004), lwd=2, lty=1, col="black")
dev.off()





# =============================================================
# Compare Three Ways of Finding c, together with Greedy Search
# =============================================================

library(rpart)
set.seed(123)
nrun <- 500
n <- 100
beta0 <- 1; beta1 <- 0; 
c <- 0; sigma <- 1
out <- matrix(0, nrun, 4)
chi2 <- matrix(0, nrun,4) ########
K0 <- 2
NLS <- F
for (i in 1:nrun){
	print(i)
	# Generate Data
	#
	# CONTINUOUS X
	# x <- runif(n, min=-1, max=1)
	# x <- rnorm(n, 0, 1)
	# 
	# CATEGORICAL X
	# x <- sample((-K0):K0, n, replace=T)/K0
	# x <- sample(c(0,1), n, replace=T)
	x <- sample(c(-1, -.5, .5, 1), n, replace=T)


	# TRANSFORM X SO THAT IT RANGES [-1, 1]
	# mid.x <- (max(x) - min(x))/2
	# x <- (x- mid.x)/(max(x)- mid.x)
	 
	y <- beta0 + beta1* sign(x<=c) + rnorm(n, mean=0, sd=sigma)
	dev0 <- (n-1)*var(y)

	# Method I: Reduction in SS - CART
	opt <- optimize(f=reduction.SS, lower =-.95, upper=.95, maximum=F, x=x, y=(y-mean(y)))
	# opt <- optimize(f=reduction.SS, lower =-.95, upper=.95, maximum=F, x=x, y=y)
	c1 <- opt$minimum 
	dev1 <-  glm(y~ sign(x <=c1), family=gaussian(link = "identity"))$deviance   # BEST!!!!!!!!!!!!!!
	# dev1 <- opt$objective 
	print(cbind(c1, dev1))
	
	# Method II: Separable LS - Profile Likelihood
	opt2 <- optimize(f=profile, lower =-.95, upper=.95, maximum=F, x=x, y=y)
	c2 <- opt2$minimum 
	dev2 <-  glm(y~ sign(x <=c2), family=gaussian(link = "identity"))$deviance 
	# dev2 <- opt2$objective


	# Method III: Nonlinear Regression  --> NOT RIGHT!!
	c3 <- 0; dev3 <- 0
	if (NLS) {
		y0 <- y - mean(y)
		a <- 100
		fit.nls <- nls(y0 ~ b1 * expit(a*(x-c)), start = list(b1 = .5, c=0), 
			control=nls.control(maxiter = 200, tol = 1e-05, minFactor = 1/1024,
            	printEval = F, warnOnly = T))
		c3 <- coef(fit.nls)[2]
		dev3 <-  glm(y~ sign(x <=c3), family=gaussian(link = "identity"))$deviance 
		# dev3 <- sum(resid(fit.nls)^2)  # Better than the option that follows  # dev3 <- (summary(fit.nls)$sigma)^2*(n-2)
	}
	
	# Greedy Search
	tree0 <- rpart(y~x, method = "anova", 
		control=rpart.control(minbucket = 10, cp = 0.000001, maxdepth = 1))
	c4 <- (tree0$splits)[4]
	z <- sign(x <= c4)
	fit.CART <- glm(y~ z, family=gaussian(link = "identity")) 
	dev4 <- fit.CART$deviance;
	
	out[i, ] <- c(c1, c2, c3, c4)
	chi2[i, ] <- dev0 - c(dev1, dev2, dev3, dev4)  # GOOD! 
	# chi2[i, ] <- c((dev0-dev1)/dev1, (dev0-dev2)/dev2, (dev0-dev3)/dev3, (dev0-dev4)/dev4)*((n-3)/2)
	# chi2[i, ] <- n*log(dev0) - n*log(c(dev1, dev2, dev3, dev4))  # WAYIII: nlog(SSE) 
}
out 
chi2

# HISTOGRAMS OF THE CUTOFF POINTS
par(mfrow=c(2,2), mar=rep(4,4))
methods <- c("Reduction SS", "Separable LS", "Nonlinear LS", "Greedy Search")
for (j in 1:ncol(out)){
	print(paste("=================", methods[j], "====================", sep=" ")) 
	cutoff <- out[,j]
	hist(cutoff, freq = FALSE, col="orangered", xlab="c", main=methods[j])
	print(cbind(mean(cutoff), sd(cutoff)))
	z <- cut(cutoff, breaks=c(-1e100, ((-4):4)/5, 1e100))
	print(table(z))
}

# DISTRIBUTION OF THE LRT 
windows()
par(mfrow=c(2,2), mar=rep(4,4))
for (j in 1:ncol(chi2)){
	print(paste("=================", methods[j], "====================", sep=" ")) 
	LRT <- chi2[,j]
	Hist(LRT, name=methods[j])
	# Hist.F(LRT, name=methods[j])
	if (j==1) legend(5, .25, legend=c(expression(chi^{2}*(1)), expression(chi^2*(2))), 
	 		col= c("red", "blue"), lty=1, lwd=2) 
}






	
# REDUCTION IN SS
expit <- function(x) {1/(1+exp(-x))}
reduction.SS <- function(c, x, y, a=200){
	if (length(x) != length(y)) stop("x and y must have the same length.")
	y <- y - mean(y)   # CENTERING Y
	# It is interesting that without centering, its performance is soooo poor. But WHY?????????  
	n <- length(x)	
	z <- expit(a*(x-c))
	# SS.YZ <- sum(y*z)
	# SS0 <- SS.YZ - sum(z*(SS.YZ/sum(z))))
	y1 <- y*z; y0 <- y*(1-z)
	SS <- var(y1)*(n-1) + var(y0)*(n-1)
	# SS1 <- sum((y1 - mean(y1))^2) + sum((y0 - mean(y0))^2)	
	# n1 <- sum(z); n0 <- n -sum(z) 	
	# ybar1 <- sum(y1)/n1; ybar0 <- sum(y0)/n0 
	# SS2 <- sum((y1 - ybar1)^2) + sum((y0 - ybar0)^2)	
	# print(rbind(SS, SS2))
	# Also, WHY THIS ABOVE METHOD DOES NOT WORK OUT???????????? 
	return(SS)
}

# PROFILE LIKELIHOOD
profile <- function(c, x, y, a=200){
	if (length(x) != length(y)) stop("x and y must have the same length.")
	y <- y - mean(y)  # Y is centered 
	z <- expit(a*(x-c))
	beta1 <- sum(z*y)/sum(z^2) 
	SS <- sum((y - beta1*z)^2) 
	return(SS)
}










#

# ======================
# EARLIER EXPERIMENTS
# ======================



n <- 100
x <- sort(runif(n))
b0 <- 1; 
b1 <- 1
sigmoid <- function(x) { 1/(1+exp(-x))}
a <- 50
y <- b0 + b1*sigmoid(a*(x-0.2))  # CART
# y <- b0 + b1*sigmoid(a*(x-0.2))*(x-0.2) # MARS
plot(x, y, type="l", lwd=2, col="green")




n <- 200
dat <- rdat(beta1=1, cutoff=.3, sd=1, n=n, discrete=F)

expit <- function(x) exp(x)/(1+exp(x))
phi <- function(c, y, x, a){
	pi <- expit(a*(x-c))
	w <- pi*(1-pi)
	r <- lsfit(x=pi, y=y)$residuals
	return(sum(r*w))
}

# PLOT OF THE ROOT-FINDING FUNCTION VS. C 
# c <- seq(0, 1, by=0.01)
c <- sort(dat$x)
phi0 <- 1:length(c)
for (i in 1:length(c)){
	phi0[i] <- phi(c[i], y=dat$y, x=dat$x, a=50)
}
plot(c, phi0, type="l", lwd=2, col="red", ylab=expression(phi(c)), xlab=expression(c)) 
abline(h=0, lwd=2, col="gray74")

# Very interesting : Uniqueness of the root is related to sample size and signal
# c <- uniroot(f=phi, y=dat$y, x=dat$x, a=50, lower = min(dat$x), upper=max(dat$x))
c <- uniroot(f=phi, y=dat$y, x=dat$x, a=50, 
	lower = quantile(dat$x, probs=.05), upper=quantile(dat$x, probs=.95))
c

# --------------------------------------------------------------------------------
# To improve this special approach, I must find all roots and quick compare their 
# corresponding SSE. 
# --------------------------------------------------------------------------------






