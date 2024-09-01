
# #################################################
# EXPERIMENT DESIGNED TO ASSESS ADAPTIVE a TO n
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
	grp <- expit(a*(x-c)); n.L <- sum(grp); sum.L <- sum(y*grp) # THE ONLY QUANTITIES THAT NEED APPROXIMATION 
	if (scale.y) SS <- sum.L^2/(n.L*(n-n.L))
	else {
		n.R <- n- n.L; sum.R <- sum(y) - sum.L
		SS <- sum.L^2/n.L + sum.R^2/n.R
	}
	return(-SS)
} 
# NOTE THAT STANDARDIZATION OF y WON'T CHANGE THE CHOICE OF c.star.

# --------------------------------------------------------------------------------------
# ALTERNATIVE OBJECTIVE FUNCTION obj.ttest() BASED ON APPROXIMATED TWO-SAMPLE t TEST 
# --------------------------------------------------------------------------------------

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




# ######################################################
# SIMULATION EXPERIMENT   - WITH a=n, sqrt(n), log(n)
# ######################################################

set.seed(888)
nrun <- 200
N <- (1:20)*50
c0 <- 0.5
OUT <- as.list(1:length(N))
# USING {rpart}
library(rpart)
minsplit <- 5
A <- c(1, 10, 30, 50, 100)
for (i in 1:length(N)) {
	n <- N[i]
	M.out <- matrix(0, nrow=nrun, ncol=length(A)+4)
	for (k in 1:nrun){	
		print(cbind(n=n, c0=c0, run=k))
		# SIMULATE DATA
		x <- runif(n, min=0, max=1)
		y <- 1 + 1*sign(x<=c0) + rnorm(n)	

		for (j in 1:length(A)){	
			# SSS WITH CONSTANT a=30 
			a <- A[j]    
			c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.01, method="ttest")
			M.out[k, j] <- c.sss 
		}

		# Greedy Search
		fit.rpart <- rpart(y~x, method="anova", control=rpart.control(minsplit = minsplit, 
			minbucket = round(minsplit /3), cp = 0, xval = 0, maxdepth = 1, 
			maxcompete = 0, maxsurrogate = 0, usesurrogate = 0))
		c.greedy <- (fit.rpart$splits)[4]
		M.out[k, length(A)+1] <- c.greedy

		# SSS WITH ADAPTIVE CHOICE OF a 
		a <- n    
		c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.01, method="ttest")
		M.out[k, length(A)+2] <- c.sss 

		# SSS WITH ADAPTIVE CHOICE OF a
		a <- sqrt(n)
		c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.01, method="ttest")
		M.out[k, length(A)+3] <- c.sss 

		# SSS WITH ADAPTIVE CHOICE OF a
		a <- log(n)
		c.sss <- bestcut.LS(x=x, y=y, a=a, scale.y=T, alpha.endcut=.01, method="ttest")
		M.out[k, length(A)+4] <- c.sss 
		}
		OUT[[i]] <- M.out
	}
}
names(OUT) <- paste("n=", N, sep="")
OUT	
save(OUT, file="OUT-an.RData")


# ---------------------------
# COMPUTE THE MSE 
# ---------------------------

c0 <- 0.5; A <- c(1, 10, 30, 50, 100)
tbl <- NULL
MSE <- function(x, c) sum((x-c)^2)/length(x)
for (i in 1:length(OUT)){
	M0 <- OUT[[i]]
	result <- apply(M0, 2, FUN=MSE, c=c0)
	tbl <- rbind(tbl, result)
}
tbl <- cbind(N, tbl)
colnames(tbl) <- c("n", paste("SSS", A, sep="."), "GS", "SSS.a=n", "SSS.sqrt.n", "SSS.ln.n")
row.names(tbl) <- NULL
tbl; dim(tbl)
tbl <- as.data.frame(tbl)

write.csv(tbl, file="Comparison-adaptive-a.csv")




# ====================
# PLOT OF MSE
# ====================


# install.packages("dplyr")
# install.packages("grid")
# install.packages("reshape2")
library(reshape2)
# install.packages("ggplot2")
# install.packages("gridExtra")
library(ggplot2)
library(gridExtra)
library(grid)
library(scales) 

tbl <- tbl.weak
colnames(tbl)
MSE <- tbl
head(MSE)
MSE0 <- data.frame(n=MSE$n, 
	# a.1 = (MSE$SSS.1-MSE$GS)/MSE$GS,
	a.10=(MSE$SSS.10-MSE$GS)/MSE$GS,
	a.30=(MSE$SSS.30-MSE$GS)/MSE$GS,
	a.50=(MSE$SSS.50-MSE$GS)/MSE$GS,
	# a.100=(MSE$SSS.100-MSE$GS)/MSE$GS,
	a.n = (MSE$"SSS.a=n"-MSE$GS)/MSE$GS,
	a.sqrt.n= (MSE$SSS.sqrt.n-MSE$GS)/MSE$GS,
	a.ln.n=(MSE$SSS.ln.n- MSE$GS)/MSE$GS)

A0 <- c(10, 30, 50, "n", "sqrt(n)", "ln(n)")
colnames(MSE0) <- c("n", paste("a = ", A0, sep=""))  

# FIGURE (a) 
# colnames(MSE)
dat0 <- MSE0
dat0 <- melt(dat0, id="n")
colnames(dat0) <- c("n", "SSS", "MSE")
fig.a <- ggplot(dat0, aes(x=n, y=MSE, color=SSS, shape=SSS)) + 
	geom_line(lwd=0.5) + geom_point() + 
	theme(legend.position="none") + 
	ggtitle(expression(paste("(a) ", beta==0.2, sep=" "))) +
	# scale_y_continuous(labels=percent) +
	labs(y = "Rel Diff in MSE vs. GS")
fig.a

tbl <- tbl.strong
colnames(tbl)
MSE <- tbl
head(MSE)
MSE0 <- data.frame(n=MSE$n, 
	# a.1 = (MSE$SSS.1-MSE$GS)/MSE$GS,
	a.10=(MSE$SSS.10-MSE$GS)/MSE$GS,
	a.30=(MSE$SSS.30-MSE$GS)/MSE$GS,
	a.50=(MSE$SSS.50-MSE$GS)/MSE$GS,
	# a.100=(MSE$SSS.100-MSE$GS)/MSE$GS,
	a.n = (MSE$"SSS.a=n"-MSE$GS)/MSE$GS,
	a.sqrt.n= (MSE$SSS.sqrt.n-MSE$GS)/MSE$GS,
	a.ln.n=(MSE$SSS.ln.n- MSE$GS)/MSE$GS)

A0 <- c(10, 30, 50, "n", expression(sqrt(n)), "ln(n)")
# A0 <- c(10, 30, 50, "n")
colnames(MSE0) <- c("n", paste("a = ", A0, sep=""), expression(a==sqrt(n)), expression(a==ln(n)))  

# FIGURE (a) 
# colnames(MSE)
dat0 <- MSE0
dat0 <- melt(dat0, id="n")
colnames(dat0) <- c("n", "SSS", "MSE")
fig.b <- ggplot(dat0, aes(x=n, y=MSE, color=SSS, shape=SSS)) + 
	geom_line(lwd=0.5) + geom_point() + 
	ggtitle(expression(paste("(b) ", beta==1.0, sep=" "))) +
	# scale_y_continuous(labels=percent) +
	labs(y = "Rel Diff in MSE vs. GS") 
# fig.b

postscript(file="fig-N06.eps", horizontal=TRUE)
multiplot(fig.a, fig.b, cols=2)
dev.off()



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#
