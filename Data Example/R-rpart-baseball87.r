

############################################################
# REGRESSION TREE EXAMPLE - THE 1987 BASEBALL SALARY DATA
############################################################

# THE 1987 BASEBALL HITTER SALARY DATA
dat <- read.table(file="bb87.dat", 
	header = F, col.names=c("id", "name", "bat86", "hit86", "hr86", "run86", 
	"rb86", "wlk86", "yrs", "batcr", "hitcr", "hrcr", "runcr", "rbcr", 
	"wlkcr", "leag86", "div86", "team86", "pos86", "puto86", "asst86", 
	"err86", "salary", "leag87", "team87", "logsalary"))
dim(dat) # 263 26
dat$leag.change <- sign(dat$leag86!=dat$leag87)
dat$team.change <- sign(dat$team86!=dat$team87)
colnames(dat)



###########################################################
# THE rpart PACKAGE USING V-FOLD CROSS-VALIDATION AND 1SE 
###########################################################

library(rpart)
colnames(dat)

set.seed(168)
control.0 <- rpart.control(minsplit=20, minbucket=5, 
	  maxdepth=15, cp=0, 
	  maxcompete=0,   						# NUMBER OF COMPETITIVE SPLITS
        maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
        xval=10)   							# SET THE VALUE V FOR V-FOLD CROSS VALIDATION
tre0 <- rpart(logsalary ~ bat86 + hit86 + hr86 + run86 + rb86 + 
	wlk86 + yrs + batcr + hitcr + hrcr + runcr + rbcr + wlkcr + 
	leag86 + div86 + team86 + pos86 + puto86 + asst86 + err86 
	+ leag87 + team87, 
	data=dat, method="anova", 
	control=control.0)


# BOOTSTRAP CART TREES TO INSPECT FOR STABILITY
if (FALSE){
	rows.b <- sample(1:nrow(dat), size=n, replace=T)
	dat.b <- dat[rows.b, ]
	dat <- dat.b
}

# CART DEFAULT TREE
# set.seed(168)
tre0 <- rpart(logsalary ~ bat86 + hit86 + hr86 + run86 + rb86 + 
	wlk86 + yrs + batcr + hitcr + hrcr + runcr + rbcr + wlkcr + 
	leag86 + div86 + team86 + pos86 + puto86 + asst86 + err86 
	+ leag87 + team87, data=dat, method="anova")

tre0; 	# THE LARGE INITIAL TREE
# summary(tre0); plot(tre0); text(tre0)

# APPLY 1-SE (1-STANDARD ERROR METHOD) FOR SELECT THE BEST TREE SIZE
# IF WE WANT TO BE LOYAL TO THE CART ENGINEERING.
printcp(tre0)
plotcp(tre0)      # 1-SE TREE SIZE SELECTION

# ----------------------------------
# OBTAINING THE BEST TREE 
# ----------------------------------

# MINIMUM CROSS-VALIDATION RELATIVE ERROR - 0SE
opt <- which.min(tre0$cptable[,"xerror"])
cp.best <- tre0$cptable[opt, "CP"]; cp.best
best.tree <- prune(tre0, cp = cp.best)
print(best.tree)
plot(best.tree, uniform=F, branch=1, compress=T, nspace=1,
     margin=.2, minbranch=.3)
text(best.tree, splits=TRUE, FUN=text, all=F,
     pretty=NULL, digits=3, use.n=T,
     fancy=FALSE)

# btre.cart.0se <- best.tree  # THE 0SE TREE GIVEN IN LOH (2002)


# ANOTHER WAY OF OBTAINING THE BEST TREE - THE 1SE APPROACH (WITH A DETAILED LOOK)
cv.error <- (tre0$cptable)[,4]
SE1 <- min(cv.error) + ((tre0$cptable)[,5])[which.min(cv.error)]      # 1SE
position <- min((1:length(cv.error))[cv.error <= SE1])
n.size  <- (tre0$cptable)[,2] + 1  # TREE SIZE IS ONE PLUS NUMBER OF SPLITS. 
best.size <- n.size[position]; best.size
best.cp <-  sqrt(tre0$cptable[position,1] *  tre0$cptable[(position-1),1])
best.cp

# SUPERIMPOSE THE BEST LAMBDA AND TREE SIZE
# -------------------------------------------
plotcp(tre0, col="red", lwd=1.5, lty=2)     # 1-SE TREE SIZE SELECTION
abline(v=best.size, lwd=2, col="green")


# ---------------------------
# OBTAIN THE OPTIMIAL SUBTREE
# ---------------------------

best.tree <- prune(tre0, cp=best.cp)
best.tree
btre.cart.1se <- best.tree
summary(best.tree)

# TREE PLOTS  AVAILABLE IN PACKAGE rpart{}
par(mfrow=c(1,1), mar=c(4, 4, 4, 4))
plot(best.tree, uniform=F, branch=1, compress=T, nspace=1,
     margin=.2, minbranch=.3)
text(best.tree, splits=TRUE, FUN=text, all=F,
     pretty=NULL, digits=3, use.n=T,
     fancy=FALSE)

# btre.cart.1se <- best.tree  # THE 1SE TREE GIVEN IN LOH (2002)



windows()
plot(best.tree, uniform=F, branch=.2, compress=F, 
     nspace=.4, margin=.4)
text(best.tree, splits=TRUE, FUN=text, all=F,
     pretty=0, digits=3, use.n=T,
     fancy=T, fwidth=0.5, fheight=0.5)   
# NOTE TAHT THE OPTION fancy=T DENOTES TERMINAL WITH A BOX  AND INTERNAL WITH ELLIPSOID. 

# POSTSRCIPT FILE FOR TREE PLOT, WHICH CAN BE FOUND IN YOUR WORKING DIRECTORY
post(best.tree, title=" ", filename = "besttree-rpart-baseball-0SE.eps", 
     pretty = TRUE, use.n = TRUE, horizontal = TRUE)

# PREDICTION AND RESIDUAL PLOTS
# -------------------------------
pred <- predict(best.tree, newdata = dat, type="vector")
resid <- dat$logsalary - pred

par(mfrow=c(2,1), mar=c(4, 6, 4, 6))
plot(dat$logsalary, pred, xlab="observed", ylab="predicted", type="p", 
	main="(a)", pch=19, col="blue", cex=.5)
abline(a=0, b=1, col="orange", lwd=2)
plot(pred, resid, xlab="predicted", ylab="residual", type="p", 
	main="(b)", col="blue", cex=.5, pch=19)
abline(h=0, lwd=2, col="orange")



################ END ##################
