#Joint model for L SNPs and K traits. (JointSum)
#Model 1: trait 1 ~ L SNPs + (K-1) traits.
#Model 2: trait 1 ~ L SNPs.

#INPUT
#B1: marginal effects on trait 1. It should be a L*1 matrix.
#S1: standard errors for marginal effects on trait 1. L*1.
#B2: marginal effects on traits to adjust for. It should be a L*(K-1) matrix.
#S2: standard errors that correspond to B2. L*(K-1).
#N: sample sizes for each coefficients in B1, B2. It should be a L*K matrix where the first column corresponds to the sample sizes for marginal effects on trait 1.
#XX: estimated covariance matrix for SNPs. Scaling this matrix will not affect the results.
#YY0: estimated correlation matrix for traits.
#adj_Y: if it is 0, adjust for SNPs only. Otherwise adjust for both SNPs and traits.
#lam: a modifying parameter in [0,1). It is used only if adj_Y=1.

#OUTPUT
#beta: coefficient estimates.
#cov: covariance matrix for coefficients.
#pvalue: p-values for coefficients.
#sigma2: estimated MSE.

JointSum=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0,adj_Y=1,lam=0){
  if(adj_Y==1){
    return(YYX4(B1,S1,B2,S2,N,XX,YY0,lam))
  }
  else{
    return(YYX5(B1,S1,B2,S2,N,XX,YY0))    
  }
}

YYX4=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0,lam=0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  
  nrx=1
  for(rx in 1:nrx){
    ny=ncol(B2)
    xx=XX[rx,rx]
    s1=S1[rx]
    b1=B1[rx]
    S21=S2[rx,]
    B21=B2[rx,]
    if(rx==1){
      yy1=N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2
    }
    else{
      yy1=yy1+N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=YY2+N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2   
    }
  }
  yy1=yy1/nrx
  YY2=YY2/nrx
  
  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix
  
  YYcross=matrix(nrow=ny+1,ncol=ny+1)
  for(i1 in 1:(ny+1)){
    for(j1 in 1:(ny+1)){
      if(i1==j1){
        YYcross[i1,j1]=YYself[i1]
      }
      else{
        YYcross[i1,j1]=YY0[i1,j1]*sqrt(YYself[i1]*YYself[j1])
      }
    }
  }
  
  nx=ncol(XX)
  A=matrix(nrow=ny+nx,ncol=ny+nx)
  A[1:nx,1:nx]=XX
  A[(nx+1):(ny+nx),1:nx]=t(XY2)
  A[1:nx,(nx+1):(ny+nx)]=XY2
  
  A[(nx+1):(ny+nx),(nx+1):(ny+nx)]=YYcross[2:(ny+1),2:(ny+1)]
  
  B=c(XY1,YYcross[1,2:(ny+1)])
  
  O=matrix(0,ncol=(ny+nx),nrow=(ny+nx))
  for(i in 1:(ny+nx)){
    O[i,i]=1
  }
  A=(A+lam*O)/(1+lam)
  
  beta=solve(A)%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx+ny))
  se=sigma2[1,1]*solve(A)   #cov
  
  if(sum(diag(se)<=0)==0){
    pvalue=c()
    for(i in 1:length(beta)){
      pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
    }
  }else{
    diag(se)=abs(diag(se))
    pvalue=beta*0-1 
  }
  
  
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}

YYX5=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  ny=ncol(B2)
  xx=XX[1,1]
  s1=S1[1]
  b1=B1[1]
  S21=S2[1,]
  B21=B2[1,]
  yy1=N[1,1]/n*(n-1)*xx*s1^2+xx*b1^2
  YY2=N[-1,1]/n*(n-1)*xx*S21^2+xx*B21^2
  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix
  nx=ncol(XX)
  A=matrix(nrow=nx,ncol=nx)
  A[1:nx,1:nx]=XX
  
  B=c(XY1)
  
  beta=solve(A)%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx))
  se=sigma2[1,1]*solve(A)   #cov
  
  pvalue=c()
  for(i in 1:length(beta)){
    pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
  }
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}