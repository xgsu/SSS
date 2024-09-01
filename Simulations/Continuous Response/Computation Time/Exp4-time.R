
rm(list=ls(all=TRUE))

# ================================================
# FUNCTIONS FOR SSS (SMOOTH SIGMOID SURROGATE)
# ================================================

# THE EXPIT FUNCTION
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED

# THE OBJECTIVE FUNCTION USED IN LEAST SQUARES WITH CONTINUOUS RESPONSE
obj0.LS <- function(c, a=50, y, x, scale.y=T){
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

# MAKE SURE Y IS SCALED SO THAT sum(y)=0
# ---------------------------------------
obj00.LS <- function(c, a=50, y, x){
	SS <- NA; n <- length(y)
	grp <- expit(a*(x-c)); n.L <- sum(grp); sum.L <- sum(y*grp)  
	sum.L^2/(n.L*(n-n.L))
}

# REDUCED SS WITHOUT CENTERING y
# --------------------------------
obj01.LS <- function(c, a=50, y, x){
	SS <- NA; n <- length(y)
	grp <- expit(a*(x-c)); n.L <- sum(grp); sum.L <- sum(y*grp)  
	n.R <- n- n.L; sum.R <- sum(y) - sum.L	
	sum.L^2/n.L + sum.R^2/n.R
} 

# THE OBJECTIVE FUNCTION USING SMOOTHED t^2 
# ------------------------------------------
obj2.LS <- function(c, a=10, y, x){
	score <- NA; n <- length(y)
	grp <- expit(a*(x-c))
	n1 <- sum(grp); n0 <- n- n1 
	y1 <- y*grp; 
	# y0 <- y*(1-grp)
	ybar1 <- sum(y1)/n1; ybar0 <- (sum(y)-sum(y1))/n0
	# print(c(n1, n0=n0, ybar1=ybar1, ybar0=ybar0))		
	# COMPUTE POOLED S2
	sp2 <- (sum(y^2) - n1*ybar1^2 - n0*ybar0^2)/(n-2)
	t <- (ybar1-ybar0)/sqrt(sp2 *(1/n1 + 1/n0)) 
	score <- t^2
	return(score)
} 




# ########################################
# Computing Time - III (with different a)
# ########################################

set.seed(123)
# SETTING II: I EXPECT TO SEE O(n^2) SPEED FOR GS IF WE MAKE K = O(n)
# WHERE K IS THE NUMBER OF DISTINCT VALUES
N <- c(2:10*10, 2:100*100)
A <- 1:100
TIME <- matrix(0, nrow=length(N), ncol=length(A)+3)
nrun <- 10
# epi <- .Machine$double.eps^0.25
epi <- 10e-3
# USING {rpart} FOR GREEDY SEARCH
library(rpart); minsplit <- 2
for (i in 1:length(N)){
	n <- N[i]; a0 <- sqrt(n)
	print(cbind(i, n)) 
	time.rpart <- time.GS <- time.SSS.sqrtn <- NULL;
	time.SSS <- matrix(0, nrow=nrun, ncol=length(A))  
	for (j in 1:nrun){
		x <- runif(n, 0, 1); mu <- 1+ 1*sign(x <= 0.5) 
		y <- mu + rnorm(n, sd=1)
		d <- .02; LB <- quantile(x, probs = d); UB <- quantile(x, probs =1-d); 	

		# METHOD I: GREEDY SEARCH USING rpart() DIRECTLY
		time.start <- proc.time() 
		# USING rpart() 
		fit <- rpart(y~x, method="anova", cp = 0, maxdepth = 1,
			control=rpart.control(minsplit = minsplit, 
			minbucket = round(minsplit/3), maxcompete = 0, maxsurrogate = 0, 
			usesurrogate = 0, xval = 0))
		time.stop <- proc.time()
		time.rpart <- c(time.rpart, sum((time.stop-time.start)[1:2]))

		# R IMPLEMETNED GREEDY SEARCH 
		time.start <- proc.time() 
		n <- length(y); wt <- rep(1, n)
		ux <- sort(unique(x));  
		y1 <- y[order(x)]; # SORT y
		temp <- cumsum(y1)[-n]
		left.wt <- cumsum(wt)[-n]
		right.wt <- sum(wt) - left.wt
		lmean <- temp/left.wt
		rmean <- -temp/right.wt
		goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y1^2)
		cuts <- (sort(x)[-n] + sort(x)[-1])/2
		c.star <- cuts[max(which.max(goodness))]
		time.stop <- proc.time()
		time.GS <- c(time.GS, sum((time.stop-time.start)[1:2]))

		# METHOD II: SSS - t2
		# a=sqrt(n)
		time.start <- proc.time() 
		cstar <- optimize(obj2.LS, lower=LB, upper=UB, maximum=T, tol=epi,
				a=sqrt(n), y=y, x=x)$maximum
		time.stop <- proc.time()
		time.SSS.sqrtn <- c(time.SSS.sqrtn, sum((time.stop-time.start)[1:2]))

		# OTHER a CHOICES
		for (k in 1:length(A)){
			time.start <- proc.time() 
			cstar <- optimize(obj2.LS, lower=LB, upper=UB, maximum=T, tol=epi,
				a=A[k], y=y, x=x)$maximum
			time.stop <- proc.time()
			time.SSS[j, k] <- sum((time.stop-time.start)[1:2])
		}
	}
	TIME[i, ] <- c(sum(time.rpart), sum(time.GS), sum(time.SSS.sqrtn), apply(time.SSS, 2, FUN=sum))
}
TIME <- as.data.frame(cbind(N, TIME))
colnames(TIME) <- c("n", "rpart", "GS", "SSS.sqrtn", A)

# TIME1 <- TIME
# TIME2 <- TIME




postscript(file="fig-N07-time.eps", width=100, height=50, horizontal=TRUE)
par(mfrow=c(1,2), mar=c(4, 5, 4,5), cex.axis=1, cex.lab=1.2, cex.main=1, cex=1)
# (a) COMPUTING TIME COMPARISON
f0 <- 1/3
n <- TIME$n; time.CART <- TIME$rpart; time.sqrtn <- TIME$"SSS.sqrtn"
plot(c(min(N), max(N)), c(0, .22), type="n", 
	xlab="n", ylab="time (seconds)", main="(a)") 
points(n, time.CART, pch=1, cex=.5, col="black") 
lines(lowess(n, time.CART, f = f0), col="black", lty=1, lwd=2)
points(n, time.sqrtn, pch=23, cex=.5, col="gray65", bg="gray65")
# lines(lowess(n, time.sqrtn, f = f0), col="gray15", lty=1, lwd=1)

rbPal <- colorRampPalette(c('red','blue'))
Col0 <- rbPal(length(A)) 
# colr <- rev(heat.colors(length(A)))
for (j in 1:length(A)){
	time.SSS <- TIME[, j+4]
	# if (j != 30) lines(lowess(n, time.SSS, f = f0), col=gray.colors(j), lty=1, lwd=0.01)  # GRAY SCALE
	lines(lowess(n, time.SSS, f = f0), col=Col0[j], lty=1, lwd=0.1)  
}
text(7000, 0.065, labels=expression(a==sqrt(n)), col="gray50", cex=1)

# (b) STEPS IN BRENT'S METHOD
load("STEP-APR2016.Rdata")
f0 <- 1/4
rbPal <- colorRampPalette(c('red','blue'))
Col0 <- rbPal(length(A)) 
# colr <- rev(heat.colors(length(A)))
plot(c(min(N), max(N)), c(7, 10.5), type="n", 
	xlab="n", ylab="steps in Brent's method", main="(b)") 
for (j in 1:length(A)){
	nstep.SSS <- STEP[, j+1]
	# lines(lowess(N, nstep.SSS, f = f0), col=gray.colors(j), lty=1, lwd=0.01)  # GRAY SCALE
	lines(lowess(N, nstep.SSS, f = f0), col=Col0[j], lty=1, lwd=0.1)  
}
text(5000, 0.18, labels=expression(a==30), col="gray50", cex=1)
legend.col(col = Col0, lev = 1:100)
mtext("Color for a Values", at=80, side=4, line=3.2, cex=1, col="black", adj=2)
dev.off()


# ADDING LEGEND WITH CONTINUOUS COLOR
legend.col <- function(col, lev){
	opar <- par
	n <- length(col)
	bx <- par("usr")
	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n
	xx <- rep(box.cx, each = 2)
	par(xpd = TRUE)
	for(i in 1:n){
		yy <- c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		polygon(xx, yy, col = col[i], border = col[i])
	}
	par(new = TRUE)
	plot(0, 0, type = "n",
	ylim = c(min(lev), max(lev)),
	yaxt = "n", ylab = "",
	xaxt = "n", xlab = "",
	frame.plot = FALSE)
	axis(side = 4, las = 2, tick = FALSE, line = .25)
	par <- opar
}



