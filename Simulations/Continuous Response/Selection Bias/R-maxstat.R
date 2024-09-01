


###############################################
# EXPERIMENT FOR ASSESSING SELECTION BIAS
###############################################

rm(list=ls(all=TRUE))
source("Functions-LS.R")

set.seed(8668)
K0 <- c(2, 3, 4, 5, 10, 20, 50, 100, 500)
p <- length(K0)
A <- c(10, 30, 50, 100)
N <- c(50, 500)
nrun <- 500
OUT <- NULL  # PREPARE DATA FOR ggplot OR LATTICE
minsplit <- 5
for (i in 1:length(N)){
	n <- N[i]
	for (j in 1:nrun){
		dat <- rdat.LS(n=n, p=length(K0), beta=c(0, 0), cutoff=.5, sd=1,
			K=K0)

		# SMOOTH SIGMOID SURROGATE (SSS)
		for (m in 1:length(A)){
			a <- A[m]
			print(cbind(n=n, run=j, a=a))
			out.tmp <- partition.LS(dat, name="1", min.ndsz=5, n0=minsplit, 
				col.y=1, col.split.var=2:(length(K0)+1), max.depth=2,  
				a=a, scale.y=T, alpha.endcut=.02, details=F)
			v.SSS <- factor.to.character((out.tmp$info)$vname)
			OUT <- rbind(OUT, c(n, j, paste("SSS a=", a, sep=""), v.SSS))  
		}
		# GREEDY SEARCH
		fit.rpart <- rpart(y~., data=dat, method="anova", 
				control=rpart.control(minsplit = minsplit, minbucket = round(minsplit /3), 
				cp = 0, xval = 1, maxdepth = 1))
		v.GS <- rownames(fit.rpart$splits)[1]
		OUT <- rbind(OUT, c(n, j, "GS", v.GS)) 
	}
}

OUT <- as.data.frame(OUT)
names(OUT) <- c("n", "run", "method", "split.var")
OUT$method <- factor(OUT$method, levels=c("SSS a=10", "SSS a=30", 
	"SSS a=50","SSS a=100", "GS"), ordered = TRUE)
OUT$n <- paste("n=", OUT$n, sep="")
OUT$n <- factor(OUT$n, levels=c("n=500", "n=50"), ordered = TRUE)
OUT$split.var <- factor(OUT$split.var, levels=paste("x", 1:p, sep=""), ordered=TRUE)
OUT <- OUT[order(OUT$n, OUT$method), ]
head(OUT)

OUT <- OUT.8668 
# library(lattice)
OUT <- OUT[OUT$method != "SSS a=100", ] 
postscript(file="fig-selectionbias-N01.eps", horizontal=T)
col0 <- c(rep("chartreuse", 3), "chocolate1", rep("chartreuse", 3), "chocolate1")
histogram(~split.var | method*n, data=OUT, layout=c(4, 2),    
	xlab=list("Splitting Variable", fontsize=16), main="", type="percent",
	ylab=list("Percentage", fontsize=16), col=col0, 
	panel=function(x, col=col, ...){
		panel.histogram(x, col=col[packet.number()], ...)
	}
) 
dev.off()



# OUT.8668 <- OUT
# OUT.123 <- OUT
# OUT.888 <- OUT





# ---------------------------------
# PLOT WITH mob RESULTS ADDED
# ---------------------------------

load("OUTmob.Rdata")

OUT <- rbind(OUT.888, OUT.mob)
OUT$method <- factor(OUT$method, levels=c("SSS a=10", "SSS a=30", 
	"SSS a=50","SSS a=100", "GS", "MOB"), ordered = TRUE)
OUT$n <- factor(OUT$n, levels=c("n=500", "n=50"), ordered = TRUE)
OUT$split.var <- factor(OUT$split.var, levels=paste("x", 1:p, sep=""), ordered=TRUE)
OUT <- OUT[order(OUT$n, OUT$method), ]


library(lattice)
postscript(file="figR-selectionbias+mob.eps", horizontal=T)
# col0 <- c(rep("darkseagreen1", 4), "goldenrod1", "cadetblue3", rep("darkseagreen1", 4), "goldenrod1", "cadetblue3")
col0 <- c(rep("chartreuse", 4), "chocolate1", "cadetblue3", rep("chartreuse", 4), "chocolate1", "cadetblue3")
histogram(~split.var | method*n, data=OUT, layout=c(6, 2),    
	xlab=list("Splitting Variable", fontsize=16), main="", type="percent",
	ylab=list("Percentage", fontsize=16), col=col0, 
	panel=function(x, col=col, ...){
		panel.histogram(x, col=col[packet.number()], ...)
	}
) 
dev.off()




# #####################################
# EXPERIMENT WITH max.t
# #####################################


# NOMINAL
n <- 500; K <- 10
y <- rnorm(n, mean=10, sd=2)
y <- scale(y)
x <- sample(1:K, size=n, replace=T)
# x <- sample(LETTERS[1:K], size=n, replace=T)  # DOESNOT WORK
dat <- data.frame(x=x, y=y)

# library(maxstat)
mod <- maxstat.test(y ~ x, data=dat, pmethod="exactGauss", minprop=0.01, maxprop=.99)
names(mod)
mod$estimate  # BEST CUTPOINT
as.numeric(mod$p.value) # P-VALUE

Nk <- tapply(rep(1,n), x, sum)
Sk <- tapply(y, x, sum)
ybars <- Sk/Nk
ord <- order(ybars)
K0 <- length(ord)
CuSum.L <- cumsum(Sk[ord])[-K0]
CuSum.R <- sum(y) - CuSum.L  
# CuSum.R <- - CuSum.L 
n.L <- cumsum(Nk[ord])[-K0]
n.R <- n - n.L
lmean <- CuSum.L/n.L
rmean <- CuSum.R/n.R
ttest.stat <-  (lmean - rmean) * sqrt(n.L*n.R/n)
max.t <- max(abs(ttest.stat))
sort(unique(x))[ord]














#

