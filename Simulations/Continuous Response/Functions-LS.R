
#####################################################
# FUNCTIONS FOR AN ALTERNATIVE TO GREEDY SEARCH
# (LS WITH CONTINUOUS DATA)
#####################################################

# install.packages("lmtest")
library(lmtest)
library(rpart)

factor.to.character <- function(fact) levels(fact)[as.numeric(fact)]


# A FUNCTION THAT SIMULATED DATA
rdat.LS <- function(n=100, p=1, beta=c(0, 2), cutoff=.5, K=rep(50,p), sd=1)
{
	X <- matrix(0, n, p)
	for (j in 1:p)  X[, j] <- sample(1:(K[j]),n, replace=T)/K[j]	
	y <- beta[1] + beta[2]*sign(X[,1]<=cutoff) + rnorm(n, mean=0, sd=sd)
	dat <- as.data.frame(cbind(y, X))
	colnames(dat) <- c("y", paste("x", 1:p, sep=""))	
   	return(dat)
}



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



# ==================================================================
# FUNCTION order.categories.IT() MAKES A CATEGORICAL VARIABLE INTO 
# ORDINAL ACCORDING TO TREATMENT EFFECTS FOR(INTERACTION TREES)
# ==================================================================

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
# y <- rnorm(100, 10, 1)
# dat <- data.frame(x=x, y=y)
# out <- order.categories.LS(dat, col.y=2, cols.cat=1, details=T)
# cbind(out$dat, x); out$OUT


# LIKELIHOOD RATIO TEST IN LINEAR REGRESSION
# USED LRT IN PACKAGE {lmtest} INSTEAD (FASTER)
LRT <- function(y, x){
	n <- length(y)
	fit1 <- lm(y~x); fit0 <- lm(y~1)
	out.tmp <- anova(fit1, fit0)
	chisq.LRT <- n*diff(log(out.tmp$RSS)) 
	chisq.LRT
}



# -------------------------------------------
# ONE SINGLE SPLIT OR PARTITION OF THE DATA
# -------------------------------------------

# WHEN USING FOR RANDOM FORESTS, SET test=NULL. 
# min.ndsz= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR CLAIMING A TERMINAL NODE
# n0= SETS THE MINIMUM NUMBER OF OBSERVATIONS FOR (n11, n10, n01, n00)
# split.var= ASKS FOR THE COLUMNS OF SPLITTING VARIABELS, INCLUDING BOTH CONTINUOUS AND CATEGORICAL ONES
# ctg= SPECIFIES THE COLUMNS OF CATEGORICAL VARIABLES
# max.depth= SPECIFIED THE MAXIMUM HEIGHT OF A TREE (ANOTHER WAY OF STOPPING THE GROWTH).
# mtry= SPECIFIES THE NUMBER OF COVARIATES IN THE RANDOMLY SELECTED SUBSETS FOR SPLITTING

partition.LS <- function(dat, test=NULL, name="1", min.ndsz=20, n0=5, 
	col.y, col.split.var, max.depth=15, mtry=length(split.var), 
	a=50, scale.y=T, alpha.endcut=.02, details=F)
{   
    # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
    out <- list(NULL) 
    name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
    n <- nrow(dat); 
    if (!is.null(test)) {n.test <- nrow(test); score.test <- NA;}  ########## TEST SAMPLE ########  
    var <- vname <- NA; cut <- NA; max.score <- -1e20;   
    y <- dat[, col.y]; vnames <- colnames(dat)
    # COMPUTE MEAN RESPONSE IN CURRENT NODE
    mean.y <- mean(y)
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
				out.tmp <- lrtest(lm(y~1), lm(y~x))
				# LRT <- out.tmp$Chisq[2]
				p.value <- out.tmp$"Pr(>Chisq)"[2]
				score <- - log10(p.value)  # logworth  
			} else {
				c.star <- bestcut.LS(x, y, a=a, scale.y=scale.y, alpha.endcut=alpha.endcut)
				grp.cstar <- expit(a*(x-c.star))
				out.tmp <- lrtest(lm(y~1), lm(y~grp.cstar))
				LRT <- out.tmp$Chisq[2]
				p.value <- pchisq(LRT, df=2, lower.tail = F)
				score <- - log10(p.value)	
			}
			if (!is.na(score) && score >= max.score) {max.score <- score; 
					var <- i; vname <- v.name; best.cut <- c.star} 
	}}}
     	if (is.null(test)) test <- dat 
	n.test <- nrow(test); y.test <- test[, col.y]; mean.y.test <- mean(y.test) 
	score <- score.test <- NA;
    	if (!is.na(var)) {
		grp <- sign(dat[,var] <= best.cut)
      	grp.test <- sign(test[,var] <= best.cut)   
		n1 <- sum(grp); n1.test <- sum(grp.test)
		if (min(n1, n-n1, n1.test, n.test-n1.test) >= n0) {
			score <- (t.test(y~grp)$statistic)^2 
			score.test <- (t.test(y.test~grp.test)$statistic)^2 
		}
	}  		
    	if (!is.na(score.test) && !is.na(score)){
        	out$name.l <- name.l; out$name.r <- name.r
        	out$left.test <- test[grp.test==1,  ]
        	out$right.test <- test[grp.test==0,  ]
            out$left  <- dat[dat[,var]<= best.cut,]
            out$right <- dat[dat[,var]> best.cut, ]
	}
    	else {var <- NA; vname <- NA; best.cut <- NA}
    	out$info <- data.frame(node=name, size=n, mean.y=mean.y,
		var = var, vname=vname, cut=best.cut, score=score, 
		score.test=score.test, size.test=n.test, mean.y.test=mean.y.test)
   	out 
}
 







