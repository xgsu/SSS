
##################################################
# FUNCTIONS IN THE R PACKAGE {exactmaxstat}
# BY 2006-09 Anne-Laure Boulesteix
##################################################



boundary<-function(x,n0,n1,c,statistic, lower=TRUE)
{
N<-n0+n1
xN<-x/N
boundary<-list()
if (statistic=="chi2")
 {
 if (lower==FALSE)
  {
  boundary$upper<-n1*x/N+n0*n1*sqrt(c)/N*sqrt(xN*(1-xN)*(1/n0+1/n1))
  return(boundary)
  }
 else
  {
  a<-n1*x/N
  b<-n0*n1*sqrt(c)/N*sqrt(xN*(1-xN)*(1/n0+1/n1))
  boundary$upper<-a+b
  boundary$lower<-a-b
  return(boundary)
  } 
 }

if (statistic=="gini")
 {
 if (lower==FALSE)
  {
  boundary$upper<-n1*x/N+x*(N-x)*sqrt(8*c/(x*(N-x)))/4
  return(boundary)
  }
 else
  {
  a<-n1*x/N
  b<-x*(N-x)*sqrt(8*c/(x*(N-x)))/4
  boundary$upper<-a+b
  boundary$lower<-a-b
  return(boundary)
  } 
 } 
 
return(boundary)
}








Fcat<-function(c,n0,n1,A,statistic)
{
N<-n0+n1
K<-length(A)
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")

number<-sum(as.numeric(permn(K,fun=sumloop_cat,n0=n0,n1=n1,A=A,c=c,statistic=statistic)))
 
return(1-number/choose(n0+n1,n1))
}

###########################
loop_cat<-function(sigma,I,k,n0,n1,A,whichboundary,c,statvector,statistic)
{
K<-length(A)
N<-n0+n1

 
upperbound<-upper(I,A,sigma)

if (k==K)
 {
 if (upperbound>=(n1-sum(I)))
  {
  number<-choose(A[K],n1-sum(I))
  }
 else
  {
  number<-0
  } 
 }
 
else 
 { 
 if (k==whichboundary)
  {
  lowerbound<-max(lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic),floor(statvector[k]+1)-sum(I))
  upperbound<-min(upperbound,n1,A[k])
  } 
 if (k<whichboundary)
  { 
  lowerbound<-max(0,lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic))
  upperbound<-min(upperbound,n1,A[k],floor(statvector[k])-sum(I))
  }
 if (k>whichboundary)
  {
  lowerbound<-max(0,lower(I,A,sigma,n0,n1,whichboundary,c,statistic=statistic))
  upperbound<-min(upperbound,n1,A[k])
  } 
 if (upperbound<lowerbound)
  {
  number<-0
  }
 else
  {  
 
number<-sum(sapply(as.list(lowerbound:upperbound),FUN=subloop_cat,I=I,k=k,n0=n0,n1=n1,A=A,sigma=sigma,whichboundary=whichboundary,c=c,statvector=statvector,statistic=statistic))
  }
 }


return(number)

}

#################
subloop_cat<-function(i,I,k,n0,n1,A,sigma,whichboundary,c,statvector,statistic)
{

number<-choose(A[k],i)*loop_cat(I=c(I,i),k=k+1,n0=n0,n1=n1,A=A,sigma=sigma,whichboundary=whichboundary,c=c,statvector=statvector,statistic=statistic)

return(number)
}


##################

upper<-function(I,A,sigma)
{
k<-length(I)+1

if (k==1)
 return(A[1])

if (sigma[k-1]>sigma[k])
 {
 upper<-ceiling(I[k-1]*A[k]/A[k-1])-1
 }

else
 {
 upper<-floor(I[k-1]*A[k]/A[k-1])
 }

return(upper)
 
}

###############

lower<-function(I,A,sigma,n0,n1,whichboundary,c,statistic)
{
k<-length(I)+1
K<-length(A)
if (k==K)
 {
 return(0)
 }
if (k>=whichboundary)
 {
 lower<-floor(A[k]*(n1-sum(I))/sum(A[k:K]))
 }
if (k<whichboundary)
 {
 n1bound<-1+floor(boundary(x=sum(A[1:whichboundary]),n0=n0,n1=n1,c=c,statistic=statistic,lower=FALSE)$upper)
 lower<-ceiling(A[k]*(n1bound-sum(I))/sum(A[k:whichboundary]))
 } 
return(lower)
}




#########################

sumloop_cat<-function(sigma,n0,n1,A,c,statistic)
{

K<-length(A)
statvector<-numeric(K-1)
sumloop_cat<-numeric(K-1)
A<-A[sigma]

for (i in 1:(K-1))
 {
 statvector[i]<-boundary(sum(A[1:i]),n0=n0,n1=n1,c=c,statistic=statistic,lower=FALSE)$upper
 }

for (k in 1:(K-1))
 {
 sumloop_cat[k]<-loop_cat(sigma,I=c(),k=1,n0=n0,n1=n1,A=A,whichboundary=k,c=c,statvector=statvector,statistic=statistic)
 }

return(sum(sumloop_cat))


}



Ford<-function(c,n0,n1,A,statistic)
{

N<-n0+n1
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")

K<-length(A)
paths<-path(c,n0,n1,A,statistic=statistic)
pathmin<-paths$pathmin
pathmax<-paths$pathmax
#A<-paths$A
Amin<-paths$Amin
Amax<-paths$Amax

xx<-c(Amin,Amax)
yy<-c(pathmin,pathmax)
x<-sort(xx)
y<-yy[order(xx)]

if (length(xx)==0)
 return(1)

b<-numeric(length(xx))

i<-2
b[1]<-choose(x[1],y[1])

while (i<=length(xx))
 {
 b[i]<-choose(x[i],y[i])-sum(sapply(1:(i-1),FUN=underpmax,i,x,y,b))
 i<-i+1
 }


pval<-0 
for (i in 1:length(xx))
 {
 pval<-pval+b[i]*choose(N-x[i],n1-y[i])
 }

pval<-1-pval/choose(N,n1)

return(pval)

}



#######################
underpmax<-function(j,i,x,y,b)
{
a<-choose(x[i]-x[j],y[i]-y[j])*b[j]
a
}

##########################



path<-function(c,n0,n1,A,statistic)
{

N<-n0+n1
K<-length(A)

pathmin<-c()
pathmax<-c()
Acum<-numeric(K-1)

quot<-1/n0+1/n1

 Amin<-c()
 Amax<-c()
 for (k in 1:(K-1))
  {
  i<-sum(A[1:k])
  Acum[k]<-i
  myboundary<-boundary(x=i,n0=n0,n1=n1,c=c,lower=TRUE,statistic=statistic)
  pmk<-floor(myboundary$upper)+1
  if (length(pathmax)>0)
   {
   a<-min(pathmax[Amax==max(Amax)])+(Acum[k]-max(Amax))-1
   }
  else
   {
   a<-n1
   }
  if (pmk<=min(i,n1,a))
   {
   pathmaxk<-seq(pmk,min(i,n1,a))
   }
  else
   {
   pathmaxk<-c()
   }

  pmk<-ceiling(myboundary$lower)-1
  if (length(pathmin>0))
   {
   a<-max(pathmin[Amin==max(Amin)])+1
   }
  else
   {
   a<-0
   }
  if (pmk>=max(0,i-n0,a))
   {
   pathmink<-seq(max(0,i-n0,a),pmk)
   }
  else
   {
   pathmink<-c()
   }

  Amin<-c(Amin,rep(Acum[k],length(pathmink)))
  Amax<-c(Amax,rep(Acum[k],length(pathmaxk)))
  pathmin<-c(pathmin,pathmink)
  pathmax<-c(pathmax,pathmaxk)
  }

list(pathmin=pathmin,pathmax=pathmax,Amin=Amin,Amax=Amax)
}


Ford2<-function(c,n0,n1,A,statistic)
{
N<-n0+n1
if (N!=sum(A))
 stop("error: you must have n0+n1=sum(A)")
K<-length(A)

if (K<=3)
 stop("when there are 3 categories, use Fcat")

b<-numeric(K-1)
b[1]<-choose(N,n1)*(1-Ford(c,n0=n0,n1=n1,A=A[c(1,K:2)],statistic=statistic))

uppervector<-numeric(N)
lowervector<-numeric(N)
for (i in 1:N)
 {
 boundaryi<-boundary(i,n0=n0,n1=n1,c=c,statistic=statistic,lower=TRUE)
 uppervector[i]<-boundaryi$upper
 lowervector[i]<-boundaryi$lower
 }

for (k in 2:(K-1))
 {
 Ak<-numeric(K)
 Ak[1:k]<-A[1:k]
 Ak[(k+1):K]<-A[K:(k+1)]

 bb<-numeric(K-k)
 # bb[i] is the number of paths crossing at k-1+i
 for (i in 1:(K-k))
  {
  bb[i]<-loop_ord2(A1=c(),k=k,n0=n0,n1=n1,A=Ak,whichboundary=k-1+i,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector)
  }
 b[k]<-sum(bb)
 }

return(1-sum(b)/choose(N,n1))
}

##################



loop_ord2<-function(A1,k,n0,n1,A,whichboundary,c,statistic,uppervector,lowervector)
{

iter<-length(A1)+1
K<-length(A)
n_p<-sum(A[as.numeric(iter>1):(iter-1)])
n1_p<-sum(A1[as.numeric(iter>1):(iter-1)])


forbidden<-c()

candidate<-max(0,n1-n1_p-sum(A[-(1:iter)])):min(A[iter],n1-n1_p)

if (iter>k)
 {
 for (i in 2:k)
  {
  if (iter>(k+1))
   {
   a<-sum(A1[c(1:(i-1),(k+1):(iter-1))])
   }
  else
   {
   a<-sum(A1[c(1:(i-1))])
   }
  u<-floor(uppervector[sum(A[c(1:(i-1),(k+1):iter)])])+1-a
  l<-ceiling(lowervector[sum(A[c(1:(i-1),(k+1):iter)])])-1-a

  if (l>=0)
   {
   forbidden<-union(forbidden,0:l)
   }
  if (u<=A[iter])
   {
   forbidden<-union(forbidden,u:A[iter])
   }
  }
 }


if (iter==K)
 {
 if (is.element(n1-n1_p,forbidden))
  {
  number<-0
  }
 else
  {
  number<-choose(A[iter],n1-n1_p)
  }
 return(number)
 }


if (iter<whichboundary)
 {
 u<-floor(uppervector[sum(A[1:iter])])+1-n1_p
 l<-ceiling(lowervector[sum(A[1:iter])])-1-n1_p
 if (l>=0)
  {
  forbidden<-union(forbidden,0:l)
  }
 if (u<=A[iter])
  {
  forbidden<-union(forbidden,u:A[iter])
  }
 }

if (iter==whichboundary)
 {
 forbidden<-union(forbidden,(ceiling(lowervector[sum(A[1:iter])])-n1_p):(floor(uppervector[sum(A[1:iter])])-n1_p))
 }



authorized<-setdiff(candidate,forbidden)
if (length(authorized)==0)
  {
  number<-0
  }
 else
  {
  number<-sum(sapply(authorized,FUN=subloop_ord2,iter=iter,A1=A1,k=k,n0=n0,n1=n1,A=A,whichboundary=whichboundary,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector))
  }
return(number)
}


#########################

subloop_ord2<-function(i,iter,A1,k,n0,n1,A,whichboundary,c,statistic,uppervector,lowervector)
{
number<-choose(A[iter],i)*loop_ord2(A1=c(A1,i),k=k,n0=n0,n1=n1,A=A,whichboundary=whichboundary,c=c,statistic=statistic,uppervector=uppervector,lowervector=lowervector)

return(number)
}





ginigain<-function(mat)
{

N<-sum(mat)
N1<-sum(mat[1,])
N2<-sum(mat[2,])
Nl<-sum(mat[,1])
Nr<-sum(mat[,2])
gain<-2/N*(N2*N1/N-mat[2,2]*mat[1,2]/Nr-mat[2,1]*mat[1,1]/Nl)
return(gain)

}





maxsel<-function(x,y=NULL,type,statistic)
{
 
statvector<-c()

 
if (!is.null(y))
 {
 notNA<-which(!is.na(x)&!is.na(y))
 x<-x[notNA]
 y<-y[notNA]
 levelsx<-sort(union(x,x))
 K<-length(levelsx)
 n0<-sum(y==0)
 n1<-length(y)-n0

 
 if (type=="ord")
  {
  for (k in 1:(K-1))
   {
   c1<-sum(y==0&(x<=levelsx[k]))
   c2<-sum(y==1&(x<=levelsx[k]))
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  }

 if (type=="ord2")
  {
  for (k in 1:(K-1))
   {
   for (j in (k+1):K)
    {
    c1<-sum(y==0&(x>levelsx[k])&(x<=levelsx[j]))
    c2<-sum(y==1&(x>levelsx[k])&(x<=levelsx[j]))
    c3<-n0-c1
    c4<-n1-c2
    if (statistic=="chi2")
     {
     statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
     }
    if (statistic=="gini")
     {
     statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
     }
    }   
   }
  }

   

 if (type=="cat")
  {
  xfac<-factor(x,labels=1:K)
  x<-as.numeric(xfac)
  prop2<-as.numeric(tapply(xfac[y==1],xfac[y==1],length))/as.numeric(tapply(xfac,xfac,length))
 
  for (k in 1:(K-1))
   {
   left<-is.element(x,order(-prop2)[1:k])
   c1<-sum(y==0&left)
   c2<-sum(y==1&left)
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  } 
 }



if (is.null(y)&ncol(as.matrix(x))>1)
 {
 if (!is.numeric(x))
  stop("x must be given as a numeric vector") 
 if (nrow(x)!=2)
  stop("x must have 2 rows and K columns")
 
 K<-ncol(x)
 zero<-c()
 for (k in 1:K)
  {
  if (sum(x[,k])==0)
   {
   zero<-c(zero,k)
   }
  }
 if (length(zero)>0) 
  {
  x<-x[,-zero]
  }
 K<-ncol(x)
 n0<-sum(x[1,])
 n1<-sum(x[2,])

 if (K<2)
  stop("x must have at least 2 columns")
 if (K<3&(type=="cat"))
  {
  type<-"ord"
  }
 
 if (type=="ord")
  {
  for (k in 1:(K-1))
   {
   c1<-sum(x[1,1:k])
   c2<-sum(x[2,1:k])
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }   
   }
  }
  

if (type=="ord2")
  {
  for (k in 1:(K-1))
   {
   for (j in (k+1):K)
    {
    c1<-sum(x[1,(k+1):j])
    c2<-sum(x[2,(k+1):j])
    c3<-n0-c1
    c4<-n1-c2
    if (statistic=="chi2")
     {
     statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
     }
    if (statistic=="gini")
     {
     statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
     }
    }   
   }
  }

  
 if (type=="cat")
  {
  prop2<-x[2,]/apply(x,FUN=sum,MARGIN=2)
  x<-x[,order(-prop2)] 
  for (k in 1:(K-1))
   {
   c1<-sum(x[1,1:k])
   c2<-sum(x[2,1:k])
   c3<-n0-c1
   c4<-n1-c2
   if (statistic=="chi2")
    {
    statvector<-c(statvector,chisq.test(matrix(c(c1,c2,c3,c4),2,2),correct=FALSE)$statistic)
    }
   if (statistic=="gini")
    {
    statvector<-c(statvector,ginigain(matrix(c(c1,c2,c3,c4),2,2)))
    }    
   }
  }  
 }
 
maxsel<-max(as.numeric(statvector))
# ALSO NEED THE CUTPOINT
return(maxsel)
}




maxsel.test<-function(x,y=NULL,type,statistic)
{

maxselcrit<-maxsel(x,y,type=type,statistic=statistic)

if (!is.null(y)&ncol(as.matrix(x))==1)
 {
 notNA<-which(!is.na(x)&!is.na(y))
 x<-x[notNA]
 y<-y[notNA]
 levelsx<-sort(union(x,x))
 K<-length(levelsx)
 n0<-sum(y==0)
 n1<-length(y)-n0
 
 A<-numeric(K)
 for (k in 1:K)
  {
  A[k]<-sum(x==levelsx[k])
  }
 }
else 
 {
 n0<-sum(x[1,])
 n1<-sum(x[2,])
 A<-apply(x,FUN=sum,MARGIN=2)
 } 

if (type=="ord")
 {
 p<-Ford(c=maxselcrit,n0=n0,n1=n1,A=A,statistic=statistic)
 return(p)
 }
 

if (type=="cat")
 {
 p<-Fcat(c=maxselcrit,n0=n0,n1=n1,A=A,statistic=statistic)
 return(p)
 }

if (type=="ord2")
 {
 p<-Ford2(c=maxselcrit,n0=n0,n1=n1,A=A,statistic=statistic)
 return(p)
 }


}



