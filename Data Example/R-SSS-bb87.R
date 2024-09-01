

############################################################
# REGRESSION TREE EXAMPLE - THE 1987 BASEBALL SALARY DATA
############################################################

rm(list=ls(all=TRUE))


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
dim(dat)

unique(sort(dat$batcr))


# ==================================================
# DATA PREPARATION: HANDLING CATEGORICAL VARIABLES
# ==================================================

source("Functions-LS.R")
cat.vnames <- c("leag86", "div86", "team86", "pos86",  "leag87", "team87")		# NAMES OF CATEGORICAL VARIABLES
cols.cat <- which(is.element(colnames(dat), cat.vnames))
out <- order.categories.LS(dat, col.y=26, cols.cat=cols.cat, details=T)
cat.info <- out$OUT
dat <- out$dat
head(dat)

# ====================================
# TEST OF FUNCTIONS 
# ====================================

col.y <- 26; 
cols.split.var <- c(3:22, 24:25, 27:28) 
vnames <- colnames(dat); vnames
vnames[cols.split.var]
mtry <- length(cols.split.var)
a <- 100; scale.y <- T; alpha.endcut <- .02; details <- F


source("Functions-LS.R")
# THE partition.LS() FUNCTION
set.seed(1288)
one.partition <- partition.LS(dat, test=NULL, name="1", min.ndsz=10, n0=5, 
	col.y, cols.split.var, max.depth=15, mtry=mtry, 
	a=50, scale.y=T, alpha.endcut=.05, details=F)
one.partition$info

# THE grow() FUNCTION - GROW A LARGE INITIAL TREE
tre0 <- grow(data=dat, test=NULL, min.ndsz=20, n0=5, 
	col.y, cols.split.var,  
	max.depth=15, mtry=mtry, 
	a=a, scale.y=scale.y, alpha.endcut=alpha.endcut, details=details)
tre0


source("Functions-LS.R")

# SPLIT-COMPLEXITY PRUNING AND BOOTSTRAP FOR BEST TREE SIZE
# set.seed(129)

# set.seed(668)
# set.seed(66188)
# DAT.SSS <- dat
dat <- DAT.SSS
boot.result <- boottrap.grow.prune(B=20, data=dat, N0=20, n0=10, 
	col.y, cols.split.var, max.depth=10, mtry=mtry, 
	a=50, global.opt.method="none", n.starts=5, 
	scale.y=T, alpha.endcut=0.02, details=F,
	LeBlanc=TRUE, min.boot.tree.size=1) 
# names(boot.result)
tree0 <- boot.result$initial.tree; tree0; 
boot.prune <- boot.result$boot.prune

n <- NROW(dat)
result <- bootstrap.size(boot.prune, tree0, n, 
	plot.it=TRUE, filename="fig-tree-size.eps", horizontal=F)
names(result)
result$G.a
result$bsize
result$btree

btre0.SSS <- result$btree$`G.log(n)`; btre0.SSS
plot.tree(btre0.SSS, textDepth=10, lines="rectangle", digits=5)

btre.SSS <- obtain.btree(tree0, bsize=5); btre.SSS 
plot.tree(btre.SSS, textDepth=6, lines="rectangle", digits=5)

# save(btre.SSS, btre.cart.1se, btre.cart.0se, file="best-trees.Rdata")

# LATEX CODES FOR THE BEST TREE STRUCTURE
plot.tree.latex(btre.SSS, file="btree-SSS-code.tex", digits=5)



# ===============================================
# UNDERSTAND WHY DIFFERENT SPLITS FROM rpart()?
# ===============================================

dat.11 <- dat[dat$batcr <= 1322, ]
t2 <- function(y, x, n0=20){
	n <- length(x)
	temp <- sort(unique(x))
	cuts <- (temp[-length(temp)] + temp[-1])/2
	t2 <- rep(0, length(cuts))
	for (i in 1:length(cuts)) {
		cut <- cuts[i]
		grp <- sign(x <= cut)
		n1 <- sum(grp); 
		if (min(n1, n-n1) >= n0) t2[i] <- (t.test(y~grp)$statistic)^2  
	}
	dat <- data.frame(x.value=cuts, t2=t2)
	dat <- dat[t2>0, ]
	return(dat)
}

par(mfrow=c(1,2), mar=rep(4,4))
t2.hitcr <- t2(y=dat.11$logsalary, x=dat.11$hitcr)
plot(t2.hitcr, type="l", col="blue", lwd=2, ylim=c(55, 170), xlab="hitcr")
abline(v=132.61006, col="orange", lwd=1.5)

t2.rbcr <- t2(y=dat.11$logsalary, x=dat.11$rbcr)
plot(t2.rbcr, type="l", col="blue", lwd=2, ylim=c(55, 170), xlab="rbcr")
abline(v=55.5, col="orange", lwd=1.5)




# ==============================================
# EXPERIMENT WITH Wilson–Hilferty approximation
# ===============================================

score <- rchisq(1000, 2)
score1 <- pmax(0, (7/9 + sqrt(2)*( (score/2)^(1/3) - 8/9))^3)
plot(score, score1, cex=.5)
hist(score, nclass=10)
mean(score)
hist(rchisq(1000,1))



# =======================================================
# 10-FOLD CROSS-VALIDATION TO COMPARE RPART AND SSS TREE
# =======================================================

# load("best-trees.Rdata")
# COMPARE "btre.cart.0se", "btre.cart.1se", AND "btre.SSS"

result$G.a
btre.SSS <- obtain.btree(tree0, bsize=5); btre.SSS 
plot.tree(btre.SSS, textDepth=6, lines="rectangle", digits=5)

dat <- DAT.SSS
dat.SSS <- send.down(dat, btre.SSS)$dat
dim(dat.SSS); head(dat.SSS)

# save(btre.cart.0se, btre.cart.1se, btre.SSS, dat.SSS, dat.bb87.rpart, file="best-trees-bb87.Rdata")
# dat.bb87.rpart["node.0se"] <- btre.cart.0se$where
# dat.bb87.rpart["node.1se"] <- btre.cart.1se$where
# head(dat.bb87.rpart)

# THE FOLLOWING ARE FAILED ATTEMPTS, WHICH WORKED WITH OLDER VERSION OF R THOUGH. 
# predict.terminal.rpart(object=btre.cart.0se, newdata=dat.bb87.rpart)
# dat.rpart.0se <- pred.rpart(btre.cart.0se, rpart.matrix(dat)) 
# rpart:::pred.rpart(btre.cart.0se, rpart:::rpart.matrix(dat))

dat <- dat.SSS;
dat["node.SSS"] <- dat$node 
dim(dat); head(dat); colnames(dat)
dat <- dat[, c(2, 26, 30)]
dat["node.0se"] <- dat.bb87.rpart["node.0se"] 
dat["node.1se"] <- dat.bb87.rpart["node.1se"] 

set.seed(286)
V <- 10
dat$cv.group <- sample(1:V, size=nrow(dat), replace =T)
dat$node <- dat$node.SSS
# dat$node <- dat$node.0se
# dat$node <- dat$node.1se
# table(dat$cv.group)
PSSE <- 0
for (v in 1:V){
	L1 <- dat[dat$cv.group!=v, ]
	L2 <- dat[dat$cv.group==v, ]
	node.mean <-  aggregate(L1$logsalary, by=list(L1$node), FUN=mean)   #### L1$node NEEDS CHANGE 
	pred.L2 <- (node.mean$x)[match(L2$node, node.mean$Group.1)]
	cbind(L2$logsalary, L2$node, pred.L2) 
	PSSE <- c(PSSE, (L2$logsalary-pred.L2)^2)
}
MPSE <- mean(PSSE); 
MPSE 

# 0.1915063
# 0.1194179
# 0.1867718


# ===============================
# DISTRIBUTION AT TERMINAL NODES
# ===============================

load("best-trees-bb87")
ls()

library(lattice)

names(dat.SSS)
names(dat.bb87.rpart)

table(dat.bb87.rpart$node.0se)
table(dat.bb87.rpart$node.1se)

h1 <- histogram(~logsalary|factor(node.0se), data=dat.bb87.rpart, 
	layout=c(8, 1), type ="density",   
	xlab="logsalary", main="(a) CART 0-SE Tree", 
	panel=function(x, ...){
		panel.histogram(x, ...)
		panel.abline(v=mean(x), col="blue", lwd=2)
		panel.densityplot(x, darg=list(bw = 0.2, kernel="gaussian"), 
			col="red", lwd=2, ...)
	}
)

h2 <- histogram(~logsalary|factor(node.1se), data=dat.bb87.rpart, 
	layout=c(6, 1), type ="density",   
	xlab="logsalary", main="(b) CART 1-SE Tree", 
	panel=function(x, ...){
		panel.histogram(x,  ...)
		panel.abline(v=mean(x), col="blue", lwd=2)
		panel.densityplot(x, darg=list(bw = 0.2, kernel="gaussian"), 
			col="red", lwd=2, ...)
	}
)

h3 <- histogram(~logsalary|node, data=dat.SSS, layout=c(5, 1), 
	type ="density",  xlab="logsalary", main="(c) SSS Tree", 
	panel=function(x, ...){
		panel.histogram(x, ...)
		panel.abline(v=mean(x), col="blue", lwd=2)
		panel.densityplot(x, darg=list(bw = 0.2, kernel="gaussian"), 
			col="red", lwd=2, ...)
	}
)


# install.packages("gridExtra")
require(gridExtra) # also loads grid

postscript(file="figR-hist-BB87.eps", horizontal=FALSE)
grid.arrange(h1, h2, h3, nrow=3, ncol=1)
dev.off()




################ END ##################
