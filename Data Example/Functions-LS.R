
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
	return(-SS)
} 
# IN FACT, STANDARDIZATION OF y WON'T CHANGE THE CHOICE OF c.star.

# SQUARED TWO-SAMPLE t TEST
obj1.LS <- function(c, a=50, y, x, scale.y=T){
	if (scale.y) y <- scale(y, center = T, scale = T)  # STANDARDIZATION OF Y MIGHT HELP WITH NUMERICAL STABILITY (PREFERRABLE)	
	score <- NA; n <- length(y)
	grp <- expit(a*(x-c))
	n1 <- sum(grp); n0 <- n-n1 
	sum.y1 <- sum(y*grp); 
	ybar1 <- sum.y1/n1; ybar0 <- (sum(y)-sum.y1)/n0	
	# COMPUTE POOLED S2
	sp2 <- (sum(y^2) - n1*ybar1^2 - n0*ybar0^2) / (n-2)
	t <- (ybar1-ybar0)/sqrt(sp2 *(1/n1 + 1/n0)) 
	score <- t^2
	return(-score)
} 


# FIND THE BEST CUTOFF POINTS FOR A CONTINUOUS variableIABLE X

# install.packages("DEoptim")
# library("DEoptim")

# GLOBAL OPTIMIZATION WILL SLOW DOWN THE PROCEDURE.
bestcut0.LS <- function(x, y, a=50, scale.y=T, alpha.endcut=.02, 
	global.opt.method = c("none", "multi.start", "global.opt"), n.starts=5)
{
	cstar <- NA
	# FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
	sigma <- sd(x); mu <- mean(x)
	x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
	LB <-  as.numeric(quantile(x, probs = alpha.endcut)); UB <- as.numeric(quantile(x, probs =1-alpha.endcut)); 
	if (LB == UB) {LB <- min(x); UB <- max(x)} 
	if (global.opt.method == "global.opt") {
		cstar <- DEoptim(obj.LS, lower=LB, upper=UB, 
			control=list(storepopfrom=1, trace=FALSE), 
			a=a, y=y, x=x, scale.y=scale.y)$optim$bestmem 	
			# print(cbind(LB, cstar, UB))
	} else if (global.opt.method == "multi.start") {
		B <- seq(LB, UB, length.out=n.starts)
		Q.min <- 1e15
		for (b in 2:n.starts) {
			OPT <- optimize(obj.LS, lower=B[b-1], upper=B[b], maximum=F, 
				a=a, y=y, x=x, scale.y=scale.y)
			if (OPT$objective < Q.min) {
				Q.min <- OPT$objective
				cstar <- OPT$minimum
			}
		}
	} else if (global.opt.method == "none") {
		cstar <- optimize(obj.LS, lower=LB, upper=UB, maximum=F, 
			a=a, y=y, x=x, scale.y=scale.y)$minimum
	} else stop("Wrong specification of gloabl.opt.method=  !!")
	cstar <- cstar*sigma + mu	# TRANSFORM BACK
	return(cstar)
} 


bestcut.LS <- function(x, y, a=50, scale.y=T, alpha.endcut=.02}
{
	cstar <- NA
	# FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
	sigma <- sd(x); mu <- mean(x)
	x <- scale(x)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
	LB <-  as.numeric(quantile(x, probs = alpha.endcut)); UB <- as.numeric(quantile(x, probs =1-alpha.endcut)); 
	if (LB == UB) {LB <- min(x); UB <- max(x)} 
	cstar <- optimize(obj.LS, lower=LB, upper=UB, maximum=F, 
			a=a, y=y, x=x, scale.y=scale.y)$minimum
	cstar <- cstar*sigma + mu	# TRANSFORM BACK
	return(cstar)
} 





# ==================================================================
# FUNCTION order.categories.IT() MAKES A CATEGORICAL variableIABLE INTO 
# ORDINAL ACCORDING TO TREATMENT EFFECTS FOR(INTERACTION TREES)
# ==================================================================

order.categories.LS <- function(dat, col.y, cols.cat, details=T){
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
		colnames(out) <- c("variableiable", "levels", "n", "mu", "rank")
		out <- as.data.frame(out)
		OUT[[j]] <- out
		x.level.ordered <- x.level[order(delta)]
		x1 <- ordered(x, levels = x.level.ordered)
		dat[, col.cat] <- as.numeric(x1)
	}
	return(list(OUT=OUT, dat=dat))
}



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
# cols.split.variables = ASKS FOR THE COLUMNS OF SPLITTING variableIABELS, INCLUDING BOTH CONTINUOUS AND CATEGORICAL ONES
# ctg= SPECIFIES THE COLUMNS OF CATEGORICAL variableIABLES
# max.depth= SPECIFIED THE MAXIMUM HEIGHT OF A TREE (ANOTHER WAY OF STOPPING THE GROWTH).
# mtry= SPECIFIES THE NUMBER OF COvariableIATES IN THE RANDOMLY SELECTED SUBSETS FOR SPLITTING

partition.LS <- function(dat, test=NULL, name="1", min.ndsz=20, n0=5, 
	col.y, cols.split.variables, max.depth=15, mtry=length(cols.split.variables), 
	a=50, global.opt.method="none", n.starts=5, 
	scale.y=T, alpha.endcut=.02, details=F)
{   
    # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL variableIABLES.
    out <- list(NULL) 
    name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
    n <- nrow(dat); 
    if (!is.null(test)) {n.test <- nrow(test); score.test <- NA;}  ########## TEST SAMPLE ########  
    variable <- vname <- NA; best.cut <- NA; max.score <- -1e20;   
    y <- dat[, col.y]; vnames <- colnames(dat)
    # COMPUTE MEAN RESPONSE IN CURRENT NODE
    mean.y <- mean(y)
    # CONTROL THE MAX TREE DEPTH
    depth <- nchar(name) 
    if (depth <= max.depth && n >= min.ndsz) {
        m.try <- ifelse(is.null(mtry), length(cols.split.variables), mtry)  
        for(i in sample(cols.split.variables, size=m.try, replace=F)) {
            x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));  
            if(length(temp) > 1) { 
			score <- NA
			if (length(temp)==2) { 
				c.star <- (temp[-length(temp)] + temp[-1])/2
				out.tmp <- lrtest(lm(y~1), lm(y~x))
				LRT <- out.tmp$Chisq[2]
				p.value <- out.tmp$"Pr(>Chisq)"[2]
				score <- - log10(p.value)  # logworth  
			} else {
				c.star <- bestcut.LS(x, y, a=a, scale.y=scale.y, alpha.endcut=alpha.endcut, 
					global.opt.method = global.opt.method, n.starts=n.starts)   ########################
				grp.cstar <- expit(a*(x-c.star))
				out.tmp <- lrtest(lm(y~1), lm(y~grp.cstar))
				LRT <- out.tmp$Chisq[2]
				p.value <- pchisq(LRT, df=2, lower.tail = F)
				score <- - log10(p.value)	
			}
			if (!is.na(score) && score >= max.score) {max.score <- LRT; 
					variable <- i; vname <- v.name; best.cut <- c.star} 
	}}}
     	if (is.null(test)) test <- dat 
	n.test <- nrow(test); y.test <- test[, col.y]; mean.y.test <- mean(y.test) 
	score <- score.test <- NA;
    	if (!is.na(variable) && !is.na(best.cut)) {
		# print(cbind(variable=variable, best.cut=best.cut))
		grp <- sign(dat[,variable] <= best.cut)
      	grp.test <- sign(test[,variable] <= best.cut)   
		n1 <- sum(grp); n1.test <- sum(grp.test)
		# print(cbind(n=n, n1=n1, n.test=n.test, n1.test=n1.test))
		if (min(n1, n-n1, n1.test, n.test-n1.test) >= n0) {
			# SQUARED t TEST
			# score <- (t.test(y~grp)$statistic)^2 
			# score.test <- (t.test(y.test~grp.test)$statistic)^2 
			# print(cbind(score, score.test))

			# LRT 
			score <- (lrtest(lm(y~1), lm(y~grp)))$Chisq[2]
			score.test <- (lrtest(lm(y.test~1), lm(y.test~grp.test)))$Chisq[2]
			# grp.cstar.test <- expit(a*(test[,variable] - best.cut))
			# print(cbind(y.test,  grp.cstar.test))
			# out.tmp <- lrtest(lm(y.test~1), lm(y.test~grp.cstar.test))
			# score.test <- out.tmp$Chisq[2]
			# print(cbind(score, score.test))

			#### TRANSFER INTO CHISQUARE 1 STATISTIC VIA Wilson–Hilferty approximation
			if (FALSE && length(unique(dat[,variable])) > 2) {
				nu <- 2
				score <- max(0, (7/9 + sqrt(nu)*((score/nu)^(1/3) - 1 + 2/(9*nu)))^3)
				score.test <- max(0, (7/9 + sqrt(nu)*((score.test/nu)^(1/3) - 1 + 2/(9*nu)))^3)
				# print(" === AFTER ====")
				# print(cbind(score, score.test))
			}
		}
	}  		
    	if (!is.na(score.test) && !is.na(score)){
        	out$name.l <- name.l; out$name.r <- name.r
        	out$left.test <- test[grp.test==1,  ]
        	out$right.test <- test[grp.test==0,  ]
            out$left  <- dat[dat[,variable]<= best.cut,]
            out$right <- dat[dat[,variable]> best.cut, ]
	}
    	else {variable <- NA; vname <- NA; best.cut <- NA}
    	out$info <- data.frame(node=name, size=n, mean.y=mean.y,
		variable = variable, vname=vname, cut=best.cut, score=score, 
		score.test=score.test, size.test=n.test, mean.y.test=mean.y.test)
   	out 
}
 


# =================================================
# THE grow.INT() FUNCTION CONSTRUCTS A LARGE TREE 
# =================================================

grow <- function(data, test=NULL, min.ndsz=20, n0=5, 
	col.y, cols.split.variables,  
	max.depth=15, mtry=length(cols.split.variables), 
	a=50, global.opt.method="none", n.starts=5, 
	scale.y=T, alpha.endcut=.02, details=F)
{
    out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
    list.nd <- list(data); 
    if (!is.null(test)) list.test <- list(test)
    name <- 1
    while (length(list.nd)!=0) {    
      for (i in 1:length(list.nd)){
	  # print(i)
        if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
		test0 <- NULL
		if (!is.null(test)) test0 <- list.test[[i]]
		split <- partition.LS(list.nd[[i]], test0, name[i], min.ndsz=min.ndsz, n0=n0, 
			col.y=col.y, cols.split.variables=cols.split.variables, max.depth=max.depth, mtry=mtry, 
			a=a, global.opt.method = global.opt.method, n.starts=n.starts, 
			scale.y=scale.y, alpha.endcut=alpha.endcut, details=details)
        	out <- rbind(out, split$info)
  	  	if (!is.null(split$left) && is.null(test)) {
			temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
			temp.name <- c(temp.name, split$name.l, split$name.r)
		
	  	} else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
			temp.list <- c(temp.list, list(split$left, split$right))
			temp.name <- c(temp.name, split$name.l, split$name.r)
			temp.test <- c(temp.test, list(split$left.test, split$right.test))
		}
	  }
	}
      list.nd <- temp.list; list.test <- temp.test; name <- temp.name
      temp.list <- temp.test <- temp.name <- NULL
    }   
    out$node <- as.character(out$node)
    out <- out[order(out$node), ]
    out
}



# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# ==========================================================


de <- function(x, tree)
{
    if(length(x) != 1) stop("The length of x in function de must be 1.")    
    y <- tree$node;  de <- NA
    if(sum(match(x, y), na.rm = T) != 0) {
        temp <- 1:length(y)
        start <- match(x, y) + 1    
        end <- length(y)
        if(start <= length(y) & nchar(y[start]) > nchar(x)) {
            temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
            if(!is.na(temp1)) end <- temp1
            de <- y[start:end]
    }}
    de
}



# ############################################################################
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ############################################################################


# =================================
# METHOD I: THE TEST SAMPLE METHOD
# =================================

# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------

prune.size.testsample <- function(tre)
{
     call <- match.call(); out <- match.call(expand = F)
     out$result <- out$size  <- out$... <- NULL
     ntest <- as.numeric(tre[1, ncol(tre)])
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre[,4])); subtree <- 1            
     a <- cbind(Ga.2=2, Ga.3=3, Ga.4=4, Ga.BIC=log(ntest))
     max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
     while (n.tmnl > 1 ) {
            # print(tre)
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
            Ga <- G - a*l 
            for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}                        
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c("variable", "vname", "cut", "score", "score.test")] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0)))     
      result <- as.data.frame(result)
      out$result <- result; out$size <- size; out$btree <- btree
      out 
}


# =================================
# METHOD II: THE BOOTSTRAP METHOD
# =================================

# -----------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE BOOTSTRAP TREE SIZE SELECTION METHOD
# -----------------------------------------------------------------------

prune.size <- function(tre)
{
     if(is.null(dim(tre))) stop("No Need to Prune Further.")
     result <- NULL; n.tmnl <- sum(is.na(tre$variable)); subtree <- 1            
     while (n.tmnl > 1 ) {
            # if (n.tmnl==5) {btre <- tre; print(btre)}
            internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            r.value <- 1:l
            for(i in 1:l) {
                branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                score <- as.numeric(as.vector(branch$score))
                r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            }
            alpha <- min(r.value)
            nod.rm <- internal[r.value == alpha]; 
            # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
            G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
            G.test <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T)
            result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                    size.tmnl=nrow(tre)-l, alpha=alpha, G=G, G.test=G.test))
            tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            tre[match(nod.rm, tre$node), c("variable", "vname", "cut", "score", "score.test")] <- NA
            n.tmnl <- sum(is.na(tre$cut))
            subtree <- subtree + 1          
      }
      # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                    size.tmnl=1, alpha=9999, G=0, G.test=0))    
      result <- as.data.frame(result)
      result
}






# ==========================================================================================  #
#  FUNCTIONS RELATED TO THE PRUNING AND THEN BOOTSTRAP FOR TREE SIZE SELECTION
# ==========================================================================================  #

# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE (TRUE) OR THE OUT-OF-BAD SAMPLE (FALSE) IN THE BOOTSTRAP PROCEDURE
# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP 

boottrap.grow.prune <- function(B=30, data, N0=20, n0=5, 
	col.y, cols.split.variables, max.depth=10, mtry=length(cols.split.variables), 
	a=50, global.opt.method="none", n.starts=5, 
	scale.y=T, alpha.endcut=0.02, details=F,
	LeBlanc=TRUE, min.boot.tree.size=1)  
{
    call <- match.call(); out <- match.call(expand = F)
    out$boot.tree <- out$boot.prune <- out$... <- NULL
    time.start <- date()
    tree0 <- grow(data=data, test=data, min.ndsz=N0, n0=n0, 
		col.y=col.y, cols.split.variables=cols.split.variables,  
		max.depth=max.depth, mtry=mtry, 
		a=a, global.opt.method=global.opt.method, n.starts=n.starts,
		scale.y=scale.y, alpha.endcut=alpha.endcut, details=details)
    print(tree0);  
    prune0 <- prune.size(tree0); 
    boot.tree <- list(tree0); boot.prune <- list(prune0) 
    b <- 1
    while (b <= B) {
        print(paste("###################### b = ", b, " ###########################", sep=""))
        # SAMPLING OBSERVATION
        samp <- sample(1:nrow(data), size=nrow(data), replace=T) 
        dat <- data[samp, ];     
        dat.oob <- data[-unique(samp),]
        n.oob <- nrow(dat.oob); # print(n.oob)        
        if (LeBlanc) {tre <- grow(data=dat, test=data, min.ndsz=N0, n0=n0, 
					col.y=col.y, cols.split.variables=cols.split.variables,  
					max.depth=max.depth, mtry=mtry, 
					a=a, global.opt.method=global.opt.method, n.starts=n.starts, 
					scale.y=scale.y, alpha.endcut=alpha.endcut, details=details)}
        else {tre <- grow(data=dat, test=dat.oob, min.ndsz=N0, n0=n0, 
					col.y=col.y, cols.split.variables=cols.split.variables,  
					max.depth=max.depth, mtry=mtry, 
					a=a, global.opt.method=global.opt.method, n.starts=n.starts,
					scale.y=scale.y, alpha.endcut=alpha.endcut, details=details)}
        print(tre)        
        if (nrow(tre)> min.boot.tree.size) {
            boot.tree <- c(boot.tree, list(tre)); 
            prune <- prune.size(tre); # print(prune)
            boot.prune <- c(boot.prune, list(prune));
            b <- b+1
        }
    }
    time.end <- date(); 
    print(paste("The Start and End time for ", B, "bootstrap runs is:"))
    print(rbind(time.start, time.end))
    out$boot.tree <- boot.tree
    out$boot.prune <- boot.prune
    # THE INITIAL LARGE TREE
    out$initial.tree <- tree0;    
    out
}   




# ========================================================================
# FUNCTION obtain.btree() OBTAINS THE BEST SUBTREE WITH KNOW SIZE bsize=
# ========================================================================

obtain.btree <- function(tre, bsize=6)
{
	if (bsize==1) { btre <- tre[1,]; btre[, c("variable", "cut", "score", "score.test")] <- NA}  
    	else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
	else {
     		n.tmnl <- sum(is.na(tre$cut)); indicator <- T            
     		while (n.tmnl > 1 && indicator ==T) {
            	if (n.tmnl==bsize) {btre <- tre; print(btre); indicator==F}
            	internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
            	r.value <- 1:l
            	for(i in 1:l) {
                		branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
                		score <- as.numeric(as.vector(branch$score))
                		r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
            	}
            	alpha <- min(r.value)
            	nod.rm <- internal[r.value == alpha]; 
            	tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
            	tre[match(nod.rm, tre$node), c("variable", "vname", "cut", "score", "score.test")] <- NA
            	n.tmnl <- sum(is.na(tre$cut))          
      	}
	}
      btre
}


# ===============================================================
# SELECT THE BEST SUBTREE SIZE AND OBTAIN THE BEST SUBTREE MODEL
# ===============================================================

bootstrap.size <- function(boot.prune, tree0, n, 
	plot.it=TRUE, filename=NULL, horizontal=T)
{   
    OUT <- as.list(NULL)
    #  COMPUTE THE ALPHA PRIME'S
    prune0 <- boot.prune[[1]] 
    n.subtree <- nrow(prune0)
    alpha <- as.numeric(as.vector(prune0$alpha));
    # temp <- c(alpha[1], alpha[-length(alpha)])  	##############
    temp <- c(0, alpha[-length(alpha)])  			############## CHANGE FIRST VALUE OF ALPHA TO 0
    alpha.prime <- sqrt(alpha*temp)  
    # cbind(alpha,  alpha.prime=prune0$alpha.prime)
    b <- length(boot.prune)
    G <- as.numeric(as.vector(prune0$G)); 
    size.tmnl <- as.numeric(as.vector(prune0$size.tmnl)); 
    subtree <- as.numeric(as.vector(prune0$subtree)); 
    # tree.penalty <- log(nrow(teeth))
    G.a <- matrix(0, n.subtree, 5)
    penalty <- c(0, 2:4, log(n))
    for (i in 1:n.subtree) {
        a.i <- alpha.prime[i]
        bias <- 0
        for (j in 2:b){
            prune.bs <- boot.prune[[j]]
            alpha.bs <- as.numeric(as.vector(prune.bs$alpha)); 
            g <- as.numeric(as.vector(prune.bs$G)); 
            g.test <- as.numeric(as.vector(prune.bs$G.test)); 
            indx <- 1
            if (sum(alpha.bs <= a.i)>0) {          
                temp1 <- which.max(which(alpha.bs<=a.i))
                indx <- ifelse(is.null(temp1), 1, temp1)
            }
            temp2 <- (g-g.test)[indx]
            bias <- bias + temp2 
            # print(cbind(i, a.i, j, bias, indx, temp2))
        }
        G.honest <- G[i] - bias/(b-1) 
        G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
    }
    out <- data.frame(cbind(size.tmnl, G.a))
    colnames(out) <- c("tmnl.size", "G", "G.2", "G.3", "G.4", "G.log(n)")
    G.a <- out

    # PLOT THE G.a WITH DIFFERENT CHOICES OF a
    if (plot.it) {
	if (!is.null(filename)) postscript(file=filename, horizontal=horizontal)
	par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS
	n.subtrees <- nrow(G.a)
	subtree.size <- G.a[,1]
	min.x <- min(subtree.size); max.x <- max(subtree.size)
	min.y <- min(G.a$"G.log(n)"); max.y <- max(G.a$G.2)
	plot(x=c(min.x, max.x), y=c(min.y, max.y), type="n", xlab="tree size", ylab="G(a)")
	for (j in 3:6) lines(subtree.size, G.a[,j], lty=j-1, col=j-1, lwd=2)
	legend(x=min.x, y=(max.y+min.y)/2, lty=2:5, col=2:5, legend=c("G(2)", "G(3)", "G(4)", "G(ln(n))")) 
	if (!is.null(filename)) dev.off()
    }   
    # OBTAIN THE BEST TREE SIZE 
    bsize <- btree <- as.list(NULL)
    Ga.cols <- c("G.2", "G.3", "G.4", "G.log(n)")
    for (j in Ga.cols) {
	best.size <- subtree.size[which.max(G.a[,j])]
	bsize[[j]] <- best.size     
      btree[[j]] <- obtain.btree(tree0, bsize=best.size)
    }
    OUT$G.a <- G.a; OUT$bsize <- bsize; OUT$btree <- btree
    return(OUT)
}	



# ==============================================================
# FUNCTION send.down() RUNS A TREE STRUCTURE DOWN A DATA SET
# ==============================================================

send.down <- function(data, tre, char.variable=1000)
{
    call <- match.call(); out <- match.call(expand = F)
    out$tree <- out$data <- out$... <- NULL
    dat <- cbind(data, node=1); tre <- cbind(tre, n.test=NA)
    cut.point <- as.vector(tre$cut); 
    split.variable <- as.numeric(as.vector(tre$variable)); 
    for (i in 1:nrow(tre)){
        in.node <- (dat$node)==(tre$node[i]);
        tre$n.test[i] <- sum(in.node)                       
        if (!is.na(split.variable[i])){
            # print(cbind(i, variable=tre$variable[i], cut=tre$cut[i]))
            variable.split <- dat[,split.variable[i]]; 
            cut <- cut.point[i]
            if (!is.element(split.variable[i], char.variable)) { 
                cut1 <- as.numeric(cut)    
                l.nd <- dat$node[in.node & variable.split <= cut1] 
                r.nd <- dat$node[in.node & variable.split > cut1]
                dat$node[in.node & variable.split <= cut1] <- paste(l.nd, 1, sep="")
                dat$node[in.node & variable.split >  cut1] <- paste(r.nd, 2, sep="")  
            }
            else {
                variable.split <- as.character(variable.split)
                cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
                l.nd <- dat$node[in.node & is.element(variable.split, cut1)] 
                r.nd <- dat$node[in.node & !is.element(variable.split, cut1)]                  
                dat$node[in.node & is.element(variable.split, cut1)] <- paste(l.nd, 1, sep="")  
                dat$node[in.node & !is.element(variable.split, cut1)] <- paste(r.nd, 2, sep="")}                   
    }}
    # print(tre)
    out$data <- dat
    out$tree <- tre
    out 
}


##############################
# TREE PLOTTING FUNCTIONS
##############################


# =============================================================
# FUNCTION plot.tree() WAS MODIFIED FROM PETER CALBOUN'S CODES
# =============================================================

# THIS FUNCTION as.numeric.factor() CONVERTS FACTOR INTO NUMERIC
as.numeric.factor <- function(x){as.numeric(levels(x))[x]}

plot.tree <- function(tree, textDepth=3, lines="rectangle", digits=5){
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes<-tree$node
  nodesBin<-gsub("1", "0", nodes)
  nodesBin<-gsub("2", "1", nodesBin)
  lastObs<-nchar(nodesBin)
  nodesBin<-substr(nodesBin,2,lastObs)
  var <- tree$vname
  cut <- as.character(tree$cut) 
  cut <- strtrim(cut, width=digits)
  size <- tree$size
  
  for(i in 1:length(nodesBin)){
    nChar<-nchar(nodesBin[i])
    if(!is.na(var[i])){
      if(lines=="rectangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      
      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1)+1/(depth+20),bquote(.(as.character(var[i]))<=.(cut[i])),cex=1)
      }
    } else {
      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1),paste("N=",size[i],sep=""),cex=1,offset=1)
      }
    }
  }
}




# =================================================================
# Generating LaTeX CODES FOR PLOTTING TREES WITH TreeTex PACKAGE 
# =================================================================

latex.plot <- function(file="tree-code.tex", tree, group=rep("I", n.node), cat.var=c(3))
{
  n.node <- nrow(tre)
  sink(file=file)
#  cat("\n \\begin{Tree} \n")
  for (i in n.node:1) {
    node <- tre[i, ] 
    if (is.na(node$variable))  {
      cat(paste("\\node{\\external\\type{circle}\\cntr{", node$size, "}\\lft{\\textbf{", group[i], "}}} \n", sep="")) 
    } else {
      if (!is.element(node$variable, cat.var)) { 
        cat(paste("\\node{\\type{frame}\\cntr{$\\texttt{V", node$variable, "} \\leq ", node$cut, "$}} \n", sep=""))
      } else {
        # cut <- node$cut
        cut1 <- unlist(strsplit(as.character(node$cut), split=" "))
        cat(paste("\\node{\\type{frame}\\cntr{\\texttt{V", node$variable, "} $ \\in \\{", paste(cut1, collapse=","), "\\} $}} \n", sep="")) 
    }}
  }
  sink()  
} 
# latex.plot(file="tree-code.tex", tree=tree0, group=rep("I", n.node), cat.var=c(3))


# =========================================================================================
# FUNCTION plot.tree.latex() GENERATES LATEX CODES FOR PLOTTING TREE IN PACKAGE pstricks
# =========================================================================================
# \usepackage{pstricks,pst-node,pst-tree}  % Plotting Trees

plot.tree.latex <- function(tree, file="tree-code.tex", digits=5)
{
	n.node <- nrow(tree)
	sink(file=file)
	# DEFINE SYMBOLS
	cat("\\begin{figure} \\centering \n")
	cat("\\newcommand{\\Tgcircle}{\\Tcircle[linecolor=black, fillcolor=gray]} \n")
	cat("\\newcommand{\\Tgoval}{\\Toval[linecolor=black, fillcolor=gray]} \n")
	cat("\\newcommand{\\Tgframe}[1]{\\Tr{\\psframebox[linecolor=black, fillstyle=solid, fillcolor=orange]{#1}}} \n")
	# OPTION
	cat("\\psset{nodesep=0.7pt, linecolor=black, treesep=1.2cm, levelsep=1.8cm} \n")
	I0 <- i0 <- NULL
	for (i in 1:n.node) {
		node.i <- as.character(tree[i, 1])
		de.i <- de(node.i, tree)
		blanks <- paste(rep(" ", (nchar(node.i)-1)*8), sep="")  # 8 SPACES IN ONE TAB		
 		n.i <- tree$size[i]
		mean.i <- ifelse(!is.null(tree$mean.y), tree$mean.y[i], " "); 	##### THIS MEASURE MAY BE DIFFERENT FOR OTHER TYPES OF TREES
		if (!is.na(de.i[1])) {	# INTERNAL NODE
			if (nchar(node.i)==1 ||  substr(node.i, nchar(node.i), nchar(node.i))=="2") 
				cat(blanks, "\\pstree{\\Tgcircle{~~}} \n", sep = "")
			else cat(blanks, "\\pstree{\\Tgcircle{~~} \\tlput{\\color{blue}", rule.i, "\\hspace{-.6in}}} \n", sep = "")				
			cat(blanks, "{ \n", sep = "") 
			I0 <- c(I0, i)
			i0 <- c(i0, i + length(de.i))
			# UPDATE THE SPLITTING RULE
			vname.i <- tree$vname[i]; cut.i <- strtrim(as.character(tree$cut[i]), width=digits); 
			rule.i <- paste("\\texttt{", vname.i, "}", "$\\leq ", cut.i, "$", sep="")
		} else if (substr(node.i, nchar(node.i), nchar(node.i))=="1") { # TERMINAL NODE
			cat(blanks, "\\Tgframe{",  mean.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{", n.i, "}}",
              		"\\tlput{\\color{blue} ", rule.i, " \\hspace{-.3in}} \n", sep = "")
		} else cat(blanks, "\\Tgframe{",  mean.i, "} \\nput{d}{\\pssucc}{\\color{black} \\textbf{",  n.i, "}} \n", sep = "")
		if (is.element(i, i0)) {
			rep0 <- rep("}", sum(i0==i))
			node.i0 <- as.character(tree[I0[i0==i][1] , 1])
			blanks0 <- paste(rep(" ", (nchar(node.i0)-1)*8), sep="")  		
			cat(blanks0, rep0, "\n", sep = "") 
		}
	}
	cat("\\end{figure} \n")
	sink()  
}
# plot.tree.latex(btre.SSS, file="tree-code.tex", digits=5)









##############################################################
##############################################################
# FUNCTION predict.terminal.rpart() PROIVDES THE TERMINAL 
# NODE MEMBERSHIP FOR EACH OBSERVATION
##############################################################
##############################################################

# IN THIS PROJECT, THIS FUNCTION IS USEFUL IN CROSS-VALIDATING CART TREE. 

library(rpart);
setOldClass(c('rpart.matrix', "matrix"))
rpart.matrix <- function(frame)
    {
    if(!inherits(frame, "data.frame"))
	    return(as.matrix(frame))
    frame$"(weights)" <- NULL
    terms <- attr(frame, "terms")
    if(is.null(terms)) predictors <- names(frame)
    else {
	a <- attributes(terms)
	predictors <- as.character(a$variables)
	removals <- NULL
	if((TT <- a$response) > 0) {
	    removals <- TT
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(TT <- a$offset)) {
	    removals <- c(removals, TT)
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(removals)) predictors <- predictors[ - removals]
        labels <- a$term.labels
	if(abs(length(labels)-length(predictors))>0)
	  predictors <- predictors[match(labels,predictors)]
	}

    factors <- sapply(frame, function(x) !is.null(levels(x)))
    characters <- sapply(frame, is.character)
    if(any(factors | characters)) {
	# change characters to factors
	for (preds in predictors[characters])
		frame[[preds]] <- as.factor(frame[[preds]])
        factors <- factors | characters
        column.levels <- lapply(frame[factors], levels)
	names(column.levels) <- (1:ncol(frame))[factors]

	# Now make them numeric
	for (preds in predictors[factors])
	     frame[[preds]] <- as.numeric(frame[[preds]])
	x <- as.matrix(frame)
	attr(x, "column.levels") <- column.levels
	}
    else x <- as.matrix(frame[predictors])
    class(x) <- "rpart.matrix"
    x
}

pred.rpart <- function(fit, x) {
    frame <- fit$frame
    nc <- frame[, c('ncompete', 'nsurrogate')]
    frame$index <- 1 + c(0, cumsum((frame$var != "<leaf>") +
                                       nc[[1]] + nc[[2]]))[-(nrow(frame)+1)]
    frame$index[frame$var == "<leaf>"] <- 0
    vnum <- match(dimnames(fit$split)[[1]], dimnames(x)[[2]])
    if (any(is.na(vnum))) stop("Tree has variables not found in new data")
    temp <- .C("pred_rpart",
		    as.integer(dim(x)),
		    as.integer(dim(frame)[1]),
		    as.integer(dim(fit$splits)),
		    as.integer(if(is.null(fit$csplit)) rep(0,2)
		               else dim(fit$csplit)),
		    as.integer(row.names(frame)),
		    as.integer(unlist(frame[,
			     c('n', 'ncompete', 'nsurrogate', 'index')])),
		    as.integer(vnum),
		    as.double(fit$splits),
		    as.integer(fit$csplit -2),
		    as.integer((fit$control)$usesurrogate),
		    as.double(x),
		    as.integer(is.na(x)),
		    where = integer(dim(x)[1]),
		    NAOK =T)
    temp <- temp$where
    names(temp) <- dimnames(x)[[1]]
    temp
 }

predict.terminal.rpart <- function(object, newdata = list()){
    if(!inherits(object, "rpart"))
	    stop("Not legitimate tree")
    if(missing(newdata))
	    where <- object$where
    else {
	if(is.null(attr(newdata, "terms"))) {
	    Terms <- delete.response(object$terms)
	    act <- (object$call)$na.action
	    if (is.null(act)) act<- na.rpart
	    newdata <- model.frame(Terms, newdata, na.action = act)
	    }
	where <- pred.rpart(object, rpart.matrix(newdata))
	}
  nodes <- row.names(object$frame)
  terminal <- nodes[where]
  return(terminal)
}









#

