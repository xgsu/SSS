
setwd("C:\\Research\\An Alternative to Greedy Search\\Programs\\Simulations\\")
rm(list=ls(all=TRUE))


library(rpart)
set.seed(123)
nrun <- 500
N <- c(50, 200)
beta0 <- 1; beta1 <- 0; 
c <- 0; sigma <- 1
out <- matrix(0, nrun, 2)
K0 <- 5
RESULT0 <- as.list(1:length(N))

for (k in 1:length(N)) {
	n <- N[k]
	for (i in 1:nrun){
		print(i)
		# Generate Data
		#
		# CONTINUOUS X
		x <- runif(n, min=-1, max=1)
		# x <- rnorm(n, 0, 1)
		# 
		# CATEGORICAL X
		# x <- sample((-K0):K0, n, replace=T)/K0

		# TRANSFORM X SO THAT IT RANGES [-1, 1]
		# mid.x <- (max(x) - min(x))/2
		# x <- (x- mid.x)/(max(x)- mid.x)
	 
		y <- beta0 + beta1* sign(x<=c) + rnorm(n, mean=0, sd=sigma)

		# Method I: Reduction in SS - CART
		opt <- optimize(f=reduction.SS, lower =-1, upper=1, maximum=F, tol =.000001, 
				x=x, y=(y-mean(y)), a=200)
		c1 <- opt$minimum 
	
		# Greedy Search
		tree0 <- rpart(y~x, method = "anova", 
			control=rpart.control(minbucket = 10, cp = 0.000001, maxdepth = 1))
		c2 <- (tree0$splits)[4]
	
		out[i, ] <- c(c1, c2)
	}
	RESULT0[[k]] <- out 
}


# CUTOFF <- as.list(1:6)
# CUTOFF[[5]] <- RESULT0[[1]]

# RESULT <- CUTOFF 

# HISTOGRAMS OF THE CUTOFF POINTS
RESULT <- RESULT0 

# HISTOGRAMS OF THE CUTOFF POINTS
par(mfrow=c(2,2), mar=rep(6,4))
methods <- c("Reduction SS", "Separable LS", "Greedy Search")
N <- c(50, 200)
breaks <- (-50:50)/50
for (k in 1:length(RESULT)){
	out <- RESULT[[k]]
	for (j in 1:ncol(out)){
		print(paste("=================", methods[j], "====================", sep=" ")) 
		cutoff <- out[,j]
		if (k==1) {hist(cutoff, probability = T, breaks=breaks, col="orangered", xlab="c", main=methods[j])}
		else {hist(cutoff, probability = T, breaks=breaks, col="orangered", xlab="c", main="")}
		if (j==ncol(out)) {mtext(paste("n=", N[k], sep=""), side=4, line =1, cex =.8)}
		print(cbind(mean(cutoff), sd(cutoff)))
		z <- cut(cutoff, breaks=c(-1e100, ((-4):4)/5, 1e100))
		print(table(z))
	}
}	



RESULT <- CUTOFF 

# postscript(file="fig2.eps", horizontal=T)
par(mfrow=c(3,4), mar=rep(4,4))
methods <- c("Smooth Sigmoid", "Greedy Search")
N <- c(50, 200)
# beta1 <- round(c(0, sqrt(3*.1), sqrt(3*1)), digit=2)
beta1 <- c("0.00", "0.55", "1.73")
breaks <- c(-1.02, (-50:50)/50, 1.02)
TBL1 <- NULL
# options(digits=4)
for (k in 1:length(RESULT)){
	out <- RESULT[[k]]
	TBL1.row <- NULL
	for (j in 1:ncol(out)){
		print(paste("=================", methods[j], "====================", sep=" ")) 
		cutoff <- out[,j]
		if (k <= 2) {
			hist(cutoff, probability = T, breaks=breaks, col="green", ylim=c(0, 20),
				xlab="c", main=methods[j])
			mtext(paste("(n=", N[j], ")", sep=""), side=3, line=0, cex=.8, col="blue")}
		else {hist(cutoff, probability = T, breaks=breaks, col="green", ylim=c(0, 20),
				xlab="c", main="")}
		if (j==ncol(out) & is.element(k, c(2,4,6))) { 
			mtext(expression(beta[1]), at=6,  side=4, line =1, cex =.8, col="blue")
			mtext(paste("=", beta1[k/2], sep=""), at = 10, side=4, line =1, cex =1, col="blue") }
		# if (k >=5) {mtext(paste("n=", N[j], sep=""), side=1, line =3, cex =1)}
		TBL1.row <- c(TBL1.row, c(mean(cutoff), sd(cutoff)))
		print(cbind(mean(cutoff), sd(cutoff)))
		z <- cut(cutoff, breaks=c(-1e100, ((-4):4)/5, 1e100))
		print(table(z))
	}
	diff.c <- out[,1] - out[,2]
	# TBL1.row <- c(beta=beta1[ceiling(k/2)], n=N[2-k%%2], TBL1.row, c(mean(diff.c), sd(diff.c)))
	TBL1.row <- c(TBL1.row, c(mean(diff.c), sd(diff.c)))
	TBL1 <- rbind(TBL1, TBL1.row)
}	
# dev.off()

# install.packages("xtable")
library(xtable)
row.names(TBL1) <- NULL
TBL1 <- as.matrix(TBL1)
x <- xtable(TBL1, caption="Summary for the Best Cutoff Points Selected by Smooth Sigmoid and Greedy Search",
	 label="tbl1", align=NULL, vsep=NULL, digits=4)
print(x, type="latex")



	
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
	# SS <- var(y1)*(n-1) + var(y0)*(n-1)
	SS <- var(y1) + var(y0)
	# SS1 <- sum((y1 - mean(y1))^2) + sum((y0 - mean(y0))^2)	
	n1 <- sum(z); n0 <- n -sum(z) 	
	# ybar1 <- sum(y1)/n1; ybar0 <- sum(y0)/n0 
	# SS2 <- sum((y1 - ybar1)^2) + sum((y0 - ybar0)^2)	
	# print(rbind(SS, SS2))
	# Also, WHY THIS ABOVE METHOD DOES NOT WORK OUT????????????
	
	# END-CUT PREFERENCE
	SS3 <- -(sum(y1))^2/n1 - (sum(y0))^2/n0 
	# SS3 <- -(sum(y1))^2 - (sum(y0))^2 
	return(SS3)
}

# ENHANCED OPTIMIZE FUNCTION

