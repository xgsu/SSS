
# #################################################
# EXPERIMENTS DESIGNED TO ASSESS SENSITIVITY TO a
# #################################################

rm(list=ls(all=TRUE))

# ================================================
# FUNCTIONS FOR SSS (SMOOTH SIGMOID SURROGATE)
# ================================================

# THE EXPIT FUNCTION
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED

# THE OBJECTIVE FUNCTION USED IN LEAST SQUARES WITH CONTINUOUS RESPONSE
obj.LS <- function(c, a=50, y, x, scale.y=T){
	SS <- NA; n <- length(y)
	grp <- expit(a*(x-c)); n.L <- sum(grp); sum.L <- sum(y*grp) # THE ONLY PLACE THAT NEEDS APPROXIMATION 
	if (scale.y) SS <- sum.L^2/(n.L*(n-n.L))
	else {
		n.R <- n- n.L; sum.R <- sum(y) - sum.L
		SS <- sum.L^2/n.L + sum.R^2/n.R
	}
	return(-SS)
} 
# IN FACT, STANDARDIZATION OF y WON'T CHANGE THE CHOICE OF c.star.


# --------------------------------------
# USING obj.ttest() SEEMS ADVANTAGEOUS
# ---------------------------------------

obj.ttest <- function(c, a=10, y, x, scale.y=T){
	if (scale.y) y <- scale(y, center = T, scale = T)  # STANDARDIZATION OF Y MIGHT HELP WITH NUMERICAL STABILITY (PREFERRABLE)
	score <- NA; n <- length(y)
	grp <- expit(a*(x-c))
	n1 <- sum(grp); n0 <- n- n1 
	y1 <- y*grp; y0 <- y*(1-grp)
	ybar1 <- sum(y1)/n1; ybar0 <- (sum(y)-sum(y1))/n0	
	sp2 <- (sum(y^2) - n1*ybar1^2 - n0*ybar0^2) / (n-2)  # COMPUTE POOLED S2
	t <- (ybar1-ybar0)/sqrt(sp2 *(1/n1 + 1/n0)) 
	score <- t^2
	return(-score)
} 


# FIND THE BEST CUTOFF POINTS FOR A CONTINUOUS VARIABLE X
bestcut.LS <- function(x, y, a=50, scale.y=T, alpha.endcut=.02, 
	method=c("ReducedSS", "ttest"), multi.start=T, n.starts=5)
{
	# FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
	sigma <- sd(x); mu <- mean(x)
	x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
	LB <- quantile(x, probs = alpha.endcut); UB <- quantile(x, probs =1-alpha.endcut); 
	if (method=="ReducedSS") obj <- obj.LS
	else obj <- obj.ttest
	if (multi.start==T) {
		B <- seq(LB, UB, length.out=n.starts)
		Q.min <- 1e15
		for (b in 2:n.starts) {
			OPT <- optimize(obj, lower=B[b-1], upper=B[b], maximum=F, 
				a=a, y=y, x=x, scale.y=scale.y)
			if (OPT$objective < Q.min) {
				Q.min <- OPT$objective
				cstar <- OPT$minimum
			}
		}
	} else {
		cstar <- optimize(obj, lower=LB, upper=UB, maximum=F, 
		a=a, y=y, x=x, scale.y=scale.y)$minimum
	}
	cstar <- cstar*sigma + mu	# TRANSFORM BACK
	return(cstar)
} 


# ##################################
# END-CUT PREFERENCE
# ##################################

set.seed(694752)
nrun <- 200
#  N <- c(50, 500)
N <- 50
A <- 1:100
OUT <- as.list(1:length(N))
# USING {rpart}
library(rpart)
minsplit <- 2
c0 <- 0.5
for (i in 1:length(N)) {
	n <- N[i]
	M.out <- matrix(0, nrow=nrun, ncol=length(A)+1)
	for (k in 1:nrun){	
		print(cbind(n=n, run=k))
		# SIMULATE DATA
		x <- runif(n, min=0, max=1)
		y <- 1 + 0*sign(x<=c0) + rnorm(n)	

		for (m in 1:length(A)){
			a <- A[m]
			# SMOOTH SIGMOID SURROGATE (SSS)
			c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.02, method="ttest")
			M.out[k, m] <- c.sss 
		}
		# Greedy Search
		fit.rpart <- rpart(y~x, method="anova", control=rpart.control(minsplit = minsplit, 
			minbucket = round(minsplit /3), cp = 0, xval = 1, maxdepth = 1))
		c.greedy <- (fit.rpart$splits)[4]
		M.out[k, length(A)+1] <- c.greedy
	}
	OUT[[i]]<- M.out
}
names(OUT) <- paste("n=", N, sep="")
OUT	

OUT.endcut <- OUT2
 

# OUT.endcut[[1]] <- OUT[[1]]
# OUT <- OUT.endcut
OUT[[2]] <- OUT.endcut[[2]]

# OUT2 <- OUT.endcut

OUT <- OUT1

postscript(file="fig03N-EndCut.eps", horizontal=T)
par(mfrow=c(1,2), mar=c(6, 4, 6, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
max.y <- c(1.25, 1.5)
c0 <- qnorm(.9)
bw0 <- 0.10
# gcol <- gray.colors(length(A), start = 0.3, end = 0.9)
# gcol <- heat.colors(length(A))
gcol <- colorRampPalette(c("blue", "red"))(length(A))
for (i in 1:length(OUT)){
	n <- ifelse(i==1, 50, 500)
	M0 <- OUT[[i]]
	d.GS <- density(M0[, ncol(M0)], bw=bw0)
	plot(d.GS, lwd=1, col="black", main=NA, xlab=expression(hat(c)), 
		xlim=c(-0.01, 1.01), ylim=c(0, max.y[i]))
	title(main=paste("(", letters[i], ") ", "n = ", n, sep=""), cex.main=1)
	polygon(d.GS, col="Ivory", border="black", lwd=.5)
	for (j in 1:(ncol(M0)-1)){
		d.SSS <- density(M0[, j], bw=bw0)
		lines(d.SSS, lwd=0.1, col=gcol[j])
	}
	if (i==1) { 
		text(0.3, 1.22, expression(a==5), col="gray40", cex=0.8)
		text(0.45, 0.82, expression(a==100), col="gray40", cex=0.8)
	}
	if (i==2) {
		text(0.38, 1.0, expression(a==5), col="gray40", cex=0.8)
		text(0.5, 0.72, expression(a==100), col="gray40", cex=0.8)
	}
}
dev.off()




# #########################################################
# WEAK/STRONG SIGNALS  - WITH VARYING a (UNIFORM X)
# #########################################################

set.seed(567)
nrun <- 200
N <- c(50, 500)
C0 <- c(.50, .80)
# C0 <- .5
A <- 1:100 
OUT <- as.list(1:(length(N)*length(C0)))
# USING {rpart}
library(rpart)
minsplit <- 2
for (i in 1:length(N)) {
	n <- N[i]
	for (j in 1:length(C0)){
		c0 <- C0[j]
		M.out <- matrix(0, nrow=nrun, ncol=length(A)+1)
		for (k in 1:nrun){	
			print(cbind(n=n, c0=c0, run=k))
			# SIMULATE DATA
			x <- runif(n, min=0, max=1)
			y <- 1 + 0.2*sign(x<=c0) + rnorm(n)	

			for (m in 1:length(A)){
				a <- A[m]
				# SMOOTH SIGMOID SURROGATE (SSS)
				c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.01, method="ttest")
				M.out[k, m] <- c.sss 
			}
			# Greedy Search
			fit.rpart <- rpart(y~x, method="anova", control=rpart.control(minsplit = minsplit , minbucket = round(minsplit /3), 
				cp = 0, xval = 1, maxdepth = 1))
			c.greedy <- (fit.rpart$splits)[4]
			M.out[k, length(A)+1] <- c.greedy
		}
		OUT[[2*(i-1)+j]] <- M.out
	}
}
names(OUT) <- c(paste("n=", 50, " c=", C0, sep=""), paste("n=", 500, " c=", C0, sep=""))
OUT	
OUT.WEAK0.2 <- OUT


# OUT.STRONG <- OUT
# OUT.WEAK0.2 <- OUT
# OUT.WEAK0.1 <- OUT




# ===============================
# DENSITY PLOT OF THE RESULTS
# ===============================

OUT <- OUT.WEAK0.2
# postscript(file="fig-N04.eps", horizontal=T)
par(mfrow=c(2,2), mar=rep(4, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
max.y <- c(2, 1.5, 2, 2)
gcol <- colorRampPalette(c("blue", "red"))(length(A))
for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
		d.GS <- density(M0[, ncol(M0)])
		plot(d.GS, lwd=1, col="black", main=NA, xlab=expression(hat(c)), 
			xlim=c(0, 1), ylim=c(0, max.y[i]))
		title(main=paste("(", letters[i], ")", sep=""), cex.main=1)
		polygon(d.GS, col="lavenderblush", border="black", lwd=.5)
		colfill <- c("black", "gray45")
		if (i==1) legend(-2, 3.0, c("GS", "SSS"), fill=colfill)
	for (j in 1:(ncol(M0)-1)){
		d.SSS <- density(M0[, j])
		lines(d.SSS, lwd=0.1, col=gcol[j])
	}
	# ADD LABEL
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	abline(v=c0, lwd=1, col="green4", lty=2)
	if (i==3) {
		# title(sub=expression(paste(c[0], "=", q[0.50],sep="")), cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", 0.5,sep="")), at=0.5, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==4) {
		# title(sub=expression(paste(c[0], "=", q[0.75],sep="")), cex.sub=1.5, col.sub="blue")
		mtext(expression(paste(c[0], "=", 0.8,sep="")), at=0.5, side=1, line=5, cex=1.5, col="blue")
	}
	if (i==2) mtext(expression(n==50), at=1.5,  side=4, line=1, cex =1.5, col="blue")
	if (i==4) mtext(expression(n==500), at=40, side=4, line=1, cex =1.5, col="blue")
}
# dev.off()


# save(OUT.endcut, OUT.STRONG, OUT.WEAK0.1, OUT.WEAK0.2, file="compasion-results.Rdata")


# ==========================
# PLOT OF ALL THE RESULTS
# ==========================

# OUT <- OUT.endcut
# OUT.STRONG 
# OUT.WEAK0.2 
# OUT.WEAK0.1

OUT <- as.list(1:10)
OUT[[1]] <- OUT.endcut[[1]]
OUT[[2]] <- OUT.endcut[[2]]
OUT[[3]] <- OUT.WEAK0.2[[1]]
OUT[[5]] <- OUT.WEAK0.2[[2]]
OUT[[4]] <- OUT.WEAK0.2[[3]]
OUT[[6]] <- OUT.WEAK0.2[[4]]
OUT[[7]] <- OUT.STRONG[[1]]
OUT[[9]] <- OUT.STRONG[[2]]
OUT[[8]] <- OUT.STRONG[[3]]
OUT[[10]] <- OUT.STRONG[[4]]


postscript(file="fig-N03.eps", horizontal=F)
N <- c(50, 500)
C0 <- c(.50, .80)
A <- 1:100 
par(mfrow=c(5,2), mar=c(2, 2.5, 2, 2.5))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
max.y <- c(1.5, 1.5, 1.5, 4, 1.6, 2.8, 8.2, 85, 7.5, 85)
# gcol <- colorRampPalette(c("green4", "Orange Red"))(length(A))
gcol <- gray.colors(length(A), start = 0.3, end = 0.9)

for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	d.GS <- density(M0[, ncol(M0)])
	plot(d.GS, lwd=1, col="black", main=NA, xlab="c", cex.axis=0.8, 
		xlim=c(0, 1), ylim=c(0.01, max.y[i]), xaxt="n")   # ifelse(i<9, "n", "s"))
	title(main=paste("(", letters[i], ")", sep=""), cex.main=0.8, line=0.4)
	polygon(d.GS, col="Old Lace", border="red", lwd=0.5)
	for (j in 1:(ncol(M0)-1)){
		d.SSS <- density(M0[, j])
		lines(d.SSS, lwd=0.1, col=gcol[length(A)-j])
	}
	# ADD LABEL
	if (i > 2){
		c0 <- ifelse(is.element(i, c(3, 4, 7, 8)), C0[1], C0[2]) 
		abline(v=c0, lwd=1, col="blue", lty=1)
	}
	if (i>=9) mtext(0:5/5, side=1, line=0, at=0:5/5, cex=0.8)
	if (i==1) mtext(expression(n==50), side=3, line=1.2, at=0.5, cex=1, col="blue")
	if (i==2) mtext(expression(n==500), side=3, line=1.2, at=0.5, cex=1, col="blue")
	if (i==2) mtext(expression(beta==0), side=4, line=1, at=0.75, cex=1, col="blue")
	if (i==4) mtext(expression(beta==0.2), side=4, line=1, at=2, cex=1, col="blue")
	if (i==6) mtext(expression(beta==0.2), side=4, line=1, at=1.4, cex=1, col="blue")
	if (i==8) mtext(expression(beta==1), side=4, line=1, at=43, cex=1, col="blue")
	if (i==10) mtext(expression(beta==1), side=4, line=1, at=43, cex=1, col="blue")
}
dev.off()


# -----------------------
# COMPUTE MSE AND PLOT
# -----------------------

tbl <-  c(paste("SSS", A, sep="."), "GS")
MSE <- function(x, c) sum((x-c)^2)/length(x)
C0 <- c(.50, .80)
for (i in 3:10){
	M0 <- OUT[[i]]
	c0 <- ifelse(is.element(i, c(3,4,7,8)), C0[1], C0[2]) 
	result <- apply(M0, 2, FUN=MSE, c=c0)
	tbl <- cbind(tbl, result)
}
tbl <- as.data.frame(tbl)
tbl
dim(tbl)


A <- 1:100; n.a <- length(A)
postscript(file="fig-N04.eps", horizontal=T)
par(mfrow=c(2,4), mar=c(2.5, 2, 2.5, 2))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
for (j in 2:ncol(tbl)){
	mse <- as.numeric(as.vector(tbl[, j]))
	y.min <- min(mse); y.max <- max(mse)
	plot(x=c(1, 100), y=c(y.min, y.max), type="n", xlab="", ylab="")
	title(main=paste("(", letters[j-1], ")", sep=""), cex.main=1, line=0.4)
	lines(A, mse[1:n.a], lty=1, col="gray45", lwd=1)
	points(A, mse[1:n.a], pch=18, col="gray45", cex=0.6)
	abline(h=mse[n.a+1], col="red", lwd=1) 
	if (j==5) mtext(expression(beta==0.2), side=4, line=1, at=0.11, cex=1.2, col="blue")
	if (j==9) mtext(expression(beta==1.0), side=4, line=1, at=0.005, cex=1.2, col="blue")
	if (j==2 || j==4) mtext(expression(n==50), side=3, line=1.2, at=50, cex=1.2, col="blue")
	if (j==3 || j==5) mtext(expression(n==500), side=3, line=1.2, at=50, cex=1.2, col="blue")
}
dev.off()








#
