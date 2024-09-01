
#####################################################
# FUNCTIONS FOR AN ALTERNATIVE TO GREEDY SEARCH
# (WITH BINARY OUTCOME)
#####################################################

# install.packages("lmtest")
# install.packages("maxstat")

source("Functions-exactmaxstat.R")


library(lmtest)
library(rpart)
library(maxstat)

# THE EXPIT FUNCTION
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED

# A FUNCTION THAT SIMULATED DATA
rdat <- function(n=100, p=1, beta=c(0, 2), cutoff=.5, K=rep(50,p))
{
	X <- matrix(0, n, p)
	for (j in 1:p)  X[, j] <- sample(1:(K[j]),n, replace=T)/K[j]	
	mu0 <- beta[1] + beta[2]*sign(X[,1]<=cutoff)
	y <- rbinom(n, size=1, prob=expit(mu0))
	dat <- as.data.frame(cbind(y, X))
	colnames(dat) <- c("y", paste("x", 1:p, sep=""))	
   	return(dat)
}



factor.to.character <- function(fact) levels(fact)[as.numeric(fact)]

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


# ==================================================================
# FUNCTION order.categories.IT() MAKES A CATEGORICAL VARIABLE INTO 
# ORDINAL ACCORDING TO TREATMENT EFFECTS FOR(INTERACTION TREES)
# ==================================================================
# THE METHOD WORKS WELL FOR BINARY RESPONSE

order.categories.LS <- function(dat, col.y, col.trt, cols.cat, details=T){
	results <- list(NULL)
	vnames <- colnames(dat)
	y <- dat[, col.y]
	p <- length(cols.cat)
	OUT <- as.list(1:p)
	names(OUT) <- colnames(dat)[cols.cat]
	for (j in 1:p){
		# j <- 1 ##############
		col.cat <- cols.cat[j]
		vname <- vnames[col.cat]
		x <- as.character(dat[, col.cat])
		x.level <- sort(unique(x))
		out <- NULL
		if (details) print(x.level)
		for (k in x.level){
			n <- sum(x==k); mu <- mean(y[x==k])
 			tmp <- c(vname, k, n, mu)
			if (details) print(tmp)
			out <- rbind(out, tmp) 
		}
		delta <- out[, ncol(out)]
		rank.level <- rank(delta, ties.method = "average")
		out <- cbind(out, rank.level)
		row.names(out) <- NULL
		colnames(out) <- c("var", "levels", "n", "mu", "rank")
		out <- as.data.frame(out)
		OUT[[j]] <- out
		x.level.ordered <- x.level[order(delta)]
		x1 <- ordered(x, levels = x.level.ordered)
		dat[, col.cat] <- as.numeric(x1)
	}
	results$OUT <- OUT
	results$dat <- dat
	return(results)
}

# EXAMPLE
# x <- rep(letters[1:4], 25)
# y <- rbinom(100, 1, .5)
# dat <- data.frame(x=x, y=y)
# out <- order.categories.LS(dat, col.y=2, cols.cat=1, details=T)
# cbind(out$dat, x); out$OUT


# -------------------------------------------
# ONE SINGLE SPLIT OR PARTITION OF THE DATA
# -------------------------------------------

# min.ndsz= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR CLAIMING A TERMINAL NODE
# n0= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR (n11, n10, n01, n00)
# col.split.var= ASKS FOR THE COLUMNS OF SPLITTING VARIABELS, INCLUDING BOTH CONTINUOUS AND CATEGORICAL ONES
# max.depth= SPECIFIED THE MAXIMUM HEIGHT OF A TREE (ANOTHER WAY OF STOPPING THE GROWTH).
# mtry= SPECIFIES THE NUMBER OF COVARIATES IN THE RANDOMLY SELECTED SUBSETS FOR SPLITTING

partition <- function(dat, name="1", min.ndsz=20, n0=5, 
	col.y, col.split.var, max.depth=15, mtry=length(split.var), 
	a=50, alpha.endcut=.02, method="entropy", multi.start=F, n.starts=5, 
	details=F)
{   
    # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
    out <- list(NULL) 
    name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
    n <- nrow(dat); 
    var <- vname <- NA; cut <- NA; max.score <- -1e20;   
    y <- dat[, col.y]; vnames <- colnames(dat)
    # COMPUTE MEAN RESPONSE IN CURRENT NODE
    p.y <- mean(y)
    # CONTROL THE MAX TREE DEPTH
    depth <- nchar(name) 
    if (depth <= max.depth && n >= min.ndsz) {
        m.try <- ifelse(is.null(mtry), length(col.split.var), mtry)  
        for(i in sample(col.split.var, size=m.try, replace=F)) {
            x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));  
            if(length(temp) > 1) { 
			score <- NA
			if (length(temp)==2) { 
				c.star <- (temp[-length(temp)] + temp[-1])/2
				out.tmp <- lrtest(glm(y~1, family=binomial(link = "logit")), 
						glm(y~x, family=binomial(link = "logit")))
				# LRT <- out.tmp$Chisq[2]
				p.value <- out.tmp$"Pr(>Chisq)"[2]
				score <- - log10(p.value)  # logworth  
			} else if (length(temp)<= 5) {
				# USE R PACKAGE exactmaxstat
				p.value <- maxsel.test(x=x, y=y, type="ord", statistic="chi2")
				########## THE exactmaxstat PACKAGE DOES NOT PUTPUT THE BEST CUTPOINT
				c.star <- mean(x)   ##################### I DON'T NEED BEST CUT FOR THIS STUDY PURPOSE THO. SO I COIN ONE FOR FUN FOR NOW. 
				score <-  - log10(p.value)
			} else {
				c.star <- bestcut.binary(x, y, a=a, alpha.endcut=alpha.endcut, 
						method=method, multi.start=multi.start, n.starts=n.starts)
				grp.cstar <- expit(a*(x-c.star))
				out.tmp <- lrtest(glm(y~1, family=binomial(link = "logit")), 
						glm(y~grp.cstar, family=binomial(link = "logit")))
				LRT <- out.tmp$Chisq[2]
				p.value <- pchisq(LRT, df=2, lower.tail = F)
				score <- - log10(p.value)	
			}
			if (!is.na(score) && score >= max.score) {max.score <- score; 
					var <- i; vname <- v.name; best.cut <- c.star} 
	}}}
    	out$info <- data.frame(node=name, size=n, p.y=p.y,
		var = var, vname=vname, cut=best.cut, score=score)
   	out 
}
 






