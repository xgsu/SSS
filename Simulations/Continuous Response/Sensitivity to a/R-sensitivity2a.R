
# #############################################################
# EXPERIMENTS DESIGNED TO ASSESS SENSITIVITY TO a - APR 2016
# #############################################################

rm(list=ls(all=TRUE))
set.seed(888)
nrun <- 100
N <- c(50, 500)
C0 <- c(.5, .8)
A <- 1:100
OUT <- as.list(1:(length(N)*length(C0)))
for (i in 1:length(N)) {
	n <- N[i]
	for (j in 1:length(C0)){
		c0 <- C0[j]
		M.out <- matrix(0, nrow=nrun, ncol=length(A))
		for (k in 1:nrun){	
			print(cbind(n=n, c0=c0, run=k))
			# SIMULATE DATA
			x <- runif(n, min=0, max=1)
			y <- 1 + 1*sign(x<=c0) + rnorm(n, mean=0, sd=1)	# SIGNAL
			
			for (l in 1:length(A)){
				# SMOOTH SIGMOID SURROGATE (SSS) 
				c.hat <- bestcut.LS(x=x, y=y, a=A[l], scale.y=T, alpha.endcut=.02)
				M.out[k, l] <- c.hat
			}
		}
		OUT[[2*(i-1)+j]] <- M.out
	}
}
names(OUT) <- c(paste("n=", 50, " c=", C0, sep=""), paste("n=", 500, " c=", C0, sep=""))
OUT	

# OUT.weak <- OUT

# ======================
# PLOT OF THE RESULTS
# ======================

postscript(file="fig-N05.eps", horizontal=TRUE)
par(mfrow=c(2,2), mar=c(6, 5, 6, 5))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
a.min <- 1

for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	mean.cut <- apply(M0, 2, mean)
	# PLOT OF CUTOFF VS. a
	plot(c(a.min, 100), c(0, 1), type="n", xlab="a", ylab=expression(hat(c)),
		main=paste("(", letters[i], ")", sep=""), cex.main=1)
	for (j in 1:nrun){
		c.hat <- M0[j,]
		lines(A, c.hat, lty=1, col="gray80", lwd=.1)
		lines(A, mean.cut, lty=1, col="green4", lwd=1)
	}
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	abline(h=c0, lwd=1, col="red", lty=2)
	if (i==3) {
		# title(sub=expression(c0==0.8), cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", 0.5,sep="")), at=50, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==4) {	
		# title(sub=expression(c0==0.8), cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", 0.8,sep="")), at=50, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==2) mtext("n=50", at=0.5,  side=4, line=1, cex =1.5, col="blue")
	if (i==4) mtext("n=500", at=0.5, side=4, line=1, cex =1.5, col="blue")
}
dev.off()







# #################################################
# EXPERIMENTS DESIGNED TO ASSESS SENSITIVITY TO a
# #################################################

rm(list=ls(all=TRUE))


set.seed(888)
nrun <- 100
N <- c(50, 500)
p0 <- c(.5, .75)
C0 <- qnorm(p0)
A <- 5:200
OUT <- as.list(1:(length(N)*length(C0)))
for (i in 1:length(N)) {
	n <- N[i]
	for (j in 1:length(C0)){
		c0 <- C0[j]
		M.out <- matrix(0, nrow=nrun, ncol=length(A))
		for (k in 1:nrun){	
			print(cbind(n=n, c0=c0, run=k))
			# SIMULATE DATA
			x <- runif(n, min=0, max=1)
			# x <- rnorm(n)
			y <- 1 + 1*sign(x<=c0) + rnorm(n, mean=0, sd=1)	# SIGNAL
			
			for (l in 1:length(A)){
				# SMOOTH SIGMOID SURROGATE (SSS) 
				c.hat <- bestcut.LS(x=x, y=y, a=A[l], scale.y=T, alpha.endcut=.02)
				M.out[k, l] <- c.hat
			}
		}
		OUT[[2*(i-1)+j]] <- M.out
	}
}
names(OUT) <- c(paste("n=", 50, " c=", C0, sep=""), paste("n=", 500, " c=", C0, sep=""))
OUT	


# ======================
# PLOT OF THE RESULTS
# ======================

postscript(file="fig3-cutoff-a-V.eps", horizontal=FALSE)
par(mfrow=c(2,2), mar=c(6, 5, 6, 5))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
a.min <- 5

for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	mean.cut <- apply(M0, 2, mean)
	# PLOT OF CUTOFF VS. a
	plot(c(a.min, 200), c(-2, 2), type="n", xlab="a", ylab=expression(hat(c)),
		main=paste("(", letters[i], ")", sep=""), cex.main=1)
	for (j in 1:nrun){
		c.hat <- M0[j,]
		lines(A, c.hat, lty=1, col="gray50", lwd=.1)
		lines(A, mean.cut, lty=3, col="gray10", lwd=1)
	}
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	abline(h=c0, lwd=1, col="red", lty=1)
	if (i==3) {
		# title(sub="c0 = 0.5", cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", q[0.50],sep="")), at=100, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==4) {	
		# title(sub="c0 = 0.75", cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", q[0.75],sep="")), at=100, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==2) mtext("n=50", at=0,  side=4, line=1, cex =1.5, col="blue")
	if (i==4) mtext("n=500", at=0, side=4, line=1, cex =1.5, col="blue")
}
dev.off()











################################################### 
# FUNCTIONS FOR SSS (SMOOTH SIGMOID SURROGATE)
###################################################



# THE EXPIT FUNCTION
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED

# THE OBJECTIVE FUNCTION USED IN LEAST SQUARES WITH CONTINUOUS RESPONSE
obj.LS <- function(c, a=50, y, x, scale.y=T){
	if (scale.y) y <- scale(y, center = T, scale = T)  # STANDARDIZATION OF Y MIGHT HELP WITH NUMERICAL STABILITY (PREFERRABLE)
	SS <- NA; n <- length(y)
	grp <- expit(a*(x-c)); n.L <- sum(grp); sum.L <- sum(y*grp) # THE ONLY PLACE THAT NEEDS APPROXIMATION 
	if (scale.y) SS <- sum.L^2/(n.L*(n-n.L))
	else {
		n.R <- n- n.L; sum.R <- sum(y) - sum.L
		SS <- sum.L^2/n.L + sum.R^2/n.R
	}
	return(SS)
} 
# IN FACT, STANDARDIZATION OF y WON'T CHANGE THE CHOICE OF c.star.


# FIND THE BEST CUTOFF POINTS FOR A CONTINUOUS VARIABLE X
bestcut.LS <- function(x, y, a=50, scale.y=T, alpha.endcut=.02){
	# FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
	sigma <- sd(x); mu <- mean(x)
	x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
	LB <- quantile(x, probs = alpha.endcut); UB <- quantile(x, probs =1-alpha.endcut); 
	cstar <- optimize(obj.LS, lower=LB, upper=UB, maximum=T, 
		a=a, y=y, x=x, scale.y=scale.y)$maximum
	cstar <- cstar*sigma + mu	# TRANSFORM BACK
	return(cstar)
} 

