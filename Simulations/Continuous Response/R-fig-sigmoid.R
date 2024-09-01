



# =====================================
# Plot of the sigmoid function
# =====================================

expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
# expit0 <- function(x, c, a=10)  (tanh(a*(x-c)/2)+1)/2 # VERIFIED
# expit <- function(x) exp(x)/(1+exp(x))
x <- -300:300/100
a <- 1:100

postscript(file="fig-sigmoid-functions.eps", horizontal=T)
par(mfrow=c(1,1), mar=c(4, 6, 4, 6))
par(cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex=1.2)
plot(c(min(x), max(x)), c(0, 1), type="n", xlab=expression(x), ylab=expression(expit(a.x)))
for (i in 1:length(a)){
	a.i <- a[i]
	y <- expit(x*a.i)
	lines(x, y, lty=1, col="gray")
}
lines(c(-3, 0, 0, 3), c(-.004, -0.004, 1.004, 1.004), lwd=2, lty=1, col="black")
text(1, 0.69, labels=c("a=1"), col="gray25")
text(1.5, 0.93, labels=c("a=2"), col="gray25")
dev.off()



# ==========================================
# Choice of a = 1, 5, 10, 20, 50, 100, 200
# ==========================================

expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
a <- c(1, 5, 10, 50, 100, 200)
x <- c(.1, .01, .001, .0001, .00001)
expit(a%*%t(x))

pnorm(0.1) - pnorm(-.1)




# ======================
# END-CUT PREFERRENCE
# ======================

expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
n <- 100
x <- 1:n/n
y <- C0 <- x
a <- 10

for (i in 1:length(C0)){
	c <- C0[i]
	# y[i] <- 1/(sum(x<=c)* sum(x>c))
	y[i] <- 1/(sum(expit(a*(x-c)))* sum(1-expit(a*(x-c))))
}

plot(C0, y, type="l", col="blue")






# ============================================================
# An Earlier Plot of the sigmoid function with a few a values
# ============================================================

expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
# expit0 <- function(x, c, a=10)  (tanh(a*(x-c)/2)+1)/2 # VERIFIED
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
