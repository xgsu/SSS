
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
obj.binary <- function(c, a=50, y, x, method=c("entropy", "Gini"))
{
	Delta.i <- NA; # REDUCTION IN IMPURITY
	n <- length(y); n1 <- sum(y==1); n0 <- n-n1
	grp <- expit(a*(x-c)); n.L <- sum(grp); n.L1 <- sum(y*grp)  #### APPROXIMATION
	n.R <- n-n.L; n.R1 <- n1 - n.L1
	
	if (method=="entropy") {
		Delta.i <- n.L1*log(n.L1) + (n.L-n.L1)*log(n.L-n.L1) + n.R1*log(n.R1) + 
			+ (n.R - n.R1)* log(n.R - n.R1) - n.L*log(n.L) - n.R*log(n.R)
	} else if (method=="Gini") {
		Delta.i <- - n.L1*(1-n.L1/n.L) - n.R1*(1-n.R1/n.R)
	}
	return(-Delta.i)
} 
# Q: do we need to set a tolerance so that min(n.L1, n.R1, n.L0, n.R0) > tolerance?? 

# FIND THE BEST CUTOFF POINTS FOR A CONTINUOUS VARIABLE X
bestcut.binary <- function(x, y, a=50, alpha.endcut=.02, 
	method=c("entropy", "Gini"), 
	multi.start=T, n.starts=5)
{
	# FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
	sigma <- sd(x); mu <- mean(x)
	x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
	LB <- quantile(x, probs = alpha.endcut); UB <- quantile(x, probs =1-alpha.endcut); 
	if (multi.start==T) {
		B <- seq(LB, UB, length.out=n.starts)
		Q.min <- 1e15
		for (b in 2:n.starts) {
			OPT <- optimize(obj.binary, lower=B[b-1], upper=B[b], maximum=F, 
				a=a, y=y, x=x, method=method)
			if (OPT$objective < Q.min) {
				Q.min <- OPT$objective
				cstar <- OPT$minimum
			}
		}
	} else {
		cstar <- optimize(obj.binary, lower=LB, upper=UB, maximum=F, 
		a=a, y=y, x=x, method=method)$minimum
	}
	cstar <- cstar*sigma + mu	# TRANSFORM BACK
	return(cstar)
} 





# ===============================
# SIMULATION EXPERIMENT
# ===============================

set.seed(88188)
nrun <- 500
N <- c(50, 500)
C0 <- c(.5, .75)
a <- 50
OUT <- as.list(1:(length(N)*length(C0)))
# USING {rpart}
library(rpart)
minsplit <- 2
for (i in 1:length(N)) {
	n <- N[i]
	for (j in 1:length(C0)){
		c0 <- C0[j]
		M.out <- matrix(0, nrow=nrun, ncol=2)
		for (k in 1:nrun){	
			print(cbind(n=n, c0=c0, run=k))
			# SIMULATE DATA
			x <- runif(n, min=0, max=1)
			# x <- rnorm(n)
			# p0 <- 1 + 1*sign(x<=c0) + 0.6*sign(x<=.2) + 0.3*sign(x<=.3)	# SIGNAL
			p0 <- 1 + 1*sign(x<=c0)	# SIGNAL
			y <- rbinom(n, size=1, prob=expit(p0))

			# SMOOTH SIGMOID SURROGATE (SSS)
			c.sss <- bestcut.binary(x=x, y=y, a=a, alpha.endcut=.02, 
				method="entropy", multi.start=F, n.starts=5)

			# Greedy Search
			fit.rpart <- rpart(y~x, method="class", parms = list(split = "information"),
				control=rpart.control(minsplit = minsplit , minbucket = round(minsplit /3), cp = 0, xval = 1, maxdepth = 1))
			c.greedy <- (fit.rpart$splits)[4]
			M.out[k, ] <- c(c.sss, c.greedy)
		}
		OUT[[2*(i-1)+j]] <- M.out
	}
}
names(OUT) <- c(paste("n=", 50, " c=", C0, sep=""), paste("n=", 500, " c=", C0, sep=""))
OUT	

# ---------------------
# COMPUTE THE MSE 
# ---------------------

tbl <-  rep(c("SSS", "greed"), 6)
MSE <- function(x, c) sum((x-c)^2)/length(x)
mean.bias <- function(x, c) mean(abs(x-c))
median.bias <- function(x, c) median(abs(x-c))

for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	result <- c(apply(M0, 2, FUN=mean),  apply(M0, 2, FUN=sd),
		apply(M0, 2, FUN=median), apply(M0, 2, FUN=mean.bias, c=c0), 
		apply(M0, 2, FUN=median.bias, c=c0), apply(M0, 2, FUN=MSE, c=c0))
	tbl <- cbind(tbl, result)
}

colnames(tbl) <- c("method", names(OUT))
row.names(tbl) <- rep(c("Mean", "SD", "Median", "Mean Bias", "Median Bias", "MSE"), rep(2,6))
tbl

write.csv(tbl, file="Comparison-Binary.csv")


# ======================
# PLOT OF THE RESULTS
# ======================

# install.packages("sm")
library(sm)

# postscript(file="fig4-SSSvsGreedy.eps", horizontal=T)
par(mfrow=c(2,2), mar=c(5, 4, 5, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	if (T){
	hist(M0[,1], prob=T, plot =T, col="darkgreen", xlim=c(0,1), 
		xlab="best cutoff", ylab="probability", 
		density=50, angle=135,
		main=paste("(", letters[i], ")", sep=""), cex.main=1)
	hist(M0[,2], prob=T, col="orange", add=T, density=10, angle=45)
	
	lines(density(M0[,1], bw = "sj"), col="darkgreen", lwd=2) 
	lines(density(M0[,2]), col="red", lwd=2) 
	}
}



# COMPARE DENSITY
# -----------------
postscript(file="fig-SSSvsGreedy-binary.eps", horizontal=T)
par(mfrow=c(2,2), mar=c(5, 4, 5, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	dat.tmp <- data.frame(cstar=c(M0[,1], M0[,2]), 
			method=factor(c(rep("sigmoid", nrun), rep("greedy", nrun))))
	sm.density.compare(dat.tmp$cstar, dat.tmp$method, xlab=expression(hat(c)), 
		ylab="density", lwd=2, lty=c(1,1), xlim=c(0, 1), alpha=.99)
	title(main=paste("(", letters[i], ")", sep=""), cex.main=1)
	colfill <- c(2:(2+length(levels(dat.tmp$method)))) 
	if (i==1) legend(.01, 2.46, levels(dat.tmp$method), fill=colfill)

	# ADD LABEL
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	abline(v=c0, lwd=1, col="gray25")
	if (i==3) title(sub="c0 = 0.5", cex.sub=1.5, line=4.2, col.sub="blue")
	if (i==4) title(sub="c0 = 0.75", cex.sub=1.5, line=4.2, col.sub="blue")
	if (i==2) mtext("n=50", at=1.0,  side=4, line=1, cex =1.5, col="blue")
	if (i==4) mtext("n=500", at=2.8, side=4, line=1, cex =1.5, col="blue")
}
dev.off()



# SMOOTH SCATTERPLOT
# ---------------------

# install.packages("KernSmooth", "RColorBrewer")
require(KernSmooth)
library(RColorBrewer)
g = 5
my.cols <- rev(brewer.pal(g, "RdYlBu"))

postscript(file="fig-SSSvsGreedy-scatter.eps", horizontal=T)
par(mfrow=c(2,2), mar=c(5, 4, 5, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
for (i in 1:length(OUT)){
	M0 <- OUT[[i]]; dim(M0)
	smoothScatter(M0[,1], M0[,2], nrpoints=.3*n, colramp=colorRampPalette(my.cols), 
		pch=19, cex=.3, col = "green1", 
		# xlab=expression(paste(hat(c), " via SSS", sep=" ")), 
		xlab="SSS", ylab="greedy") 
	title(main=paste("(", letters[i], ")", sep=""), cex.main=1)

	# ADD LABEL
	c0 <- ifelse(is.element(i, c(1, 3)), C0[1], C0[2]) 
	abline(v=c0, lwd=1, col="gray95")
	abline(h=c0, lwd=1, col="gray95")
	if (i==3) title(sub="c0 = 0.5", cex.sub=1.5, line=4.2, col.sub="blue")
	if (i==4) title(sub="c0 = 0.75", cex.sub=1.5, line=4.2, col.sub="blue")
	if (i==2) mtext("n=50", at=.5,  side=4, line=1, cex =1.5, col="blue")
	if (i==4) mtext("n=500", at=.5, side=4, line=1, cex =1.5, col="blue")
}
dev.off()






# ==============================
# PLOT USED IN THE ARTICLE
# ==============================

require(KernSmooth)
library(RColorBrewer)
g = 5
my.cols <- rev(brewer.pal(g, "RdYlBu"))

postscript(file="fig-comparison-binary.eps", horizontal=T)
par(mfrow=c(2,2), mar=c(5, 4, 5, 4))
par(cex.axis=1, cex.lab=1, cex.main=1, cex=1)
for (i in c(1, 3)){
	M0 <- OUT[[i]]; 
	smoothScatter(M0[,1], M0[,2], nrpoints=.3*n, colramp=colorRampPalette(my.cols), 
		pch=19, cex=.3, col = "green1", 
		# xlab=expression(paste(hat(c), " via SSS", sep=" ")), 
		xlab="sigmoid", ylab="greedy") 
	title(main=paste("(", letters[i], ")", sep=""), cex.main=1)
	abline(v=0.5, lwd=1, col="gray95")
	abline(h=0.5, lwd=1, col="gray95")


	dat.tmp <- data.frame(cstar=c(M0[,1], M0[,2]), 
			method=factor(c(rep("sigmoid", nrun), rep("greedy", nrun))))
	sm.density.compare(dat.tmp$cstar, dat.tmp$method, xlab=expression(hat(c)), 
		ylab="density", lwd=2, lty=c(1,1), xlim=c(0, 1), alpha=.99)
	title(main=paste("(", letters[i+1], ")", sep=""), cex.main=1)
	colfill <- c(2:(2+length(levels(dat.tmp$method)))) 
	if (i==1) legend(.01, 2.46, levels(dat.tmp$method), fill=colfill)
	abline(v=0.5, lwd=.6, col="gray25")

	# ADD LABEL 

	if (i==1) mtext("n=50", at=1.25,  side=4, line=1, cex =1.5, col="blue")
	if (i==3) mtext("n=500", at=4, side=4, line=1, cex =1.5, col="blue")
}
dev.off()


