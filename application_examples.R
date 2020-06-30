source('functions.R')

require(rlist) 
require(MASS)
require(mvtnorm)

n=500
p<-4
m<-2

X<-array(0,c(p,m,n))   
for(i in 1:n) {
		aux<-matrix(1,2,4)
    aux[,2]<-runif(2,0,10)
    aux[,3]<-runif(2,0,1)
    aux[,4]<-runif(2,0,5)
    X[,,i]<-t(aux)
}

max.iter=500
prec=1e-3

# True values
betas<-matrix(c(5,1,0.75,-0.20),ncol=1)
Sigmas<-matrix(c(0.5,0,0,2),ncol=2)
lambdas<-matrix(c(-3,3),ncol=1)
nu=5   # Slash and T


y=matrix(0,n,m)
for(i in 1:n) {
     mu<-t(X[,,i])%*%betas
     y[i,]<-rsmsnc(n=1,mu,Sigmas,lambdas,dist="STEC",nu=nu)
}

# SNC     
thetaSNC=EM.SNC(y,X,max.iter=max.iter,prec=prec) 
beta0=thetaSNC$beta
Sigma0=thetaSNC$Sigma
eta0=thetaSNC$eta
   
# STEC
thetaSTEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="STEC",max.iter=max.iter,prec=prec) 

# SSLEC
thetaSSLEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="SSLEC",max.iter=max.iter,prec=prec) 
        
# SCEC
thetaSCEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="SCEC",max.iter=max.iter,prec=prec) 


#############################################################################################################

# Application - real data set

dados=read.table("australia3.txt",header=T)
attach(dados)
names(dados)
                                                                                     
n=nrow(dados)
p=2

y=cbind(Hg,ssf)

q0=3
m=2*q0


X<-array(0,c(m,p,n))
for(i in 1:n){
    auxi1=cbind(1,Hc[i],Fe[i])
    auxi2=cbind(1,Bfat[i],lbm[i])
    mi=matrix(0,m,p)
    mi[1:q0,1]=auxi1
    mi[(q0+1):m,2]=auxi2
    X[,,i]=mi  
}
    
# SNC     
thetaSNC=EM.SNC(y,X,max.iter=max.iter,prec=prec) 
beta0=thetaSNC$beta
Sigma0=thetaSNC$Sigma
eta0=thetaSNC$eta
   
# STEC
thetaSTEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="STEC",max.iter=max.iter,prec=prec) 

# SSLEC
thetaSSLEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="SSLEC",max.iter=max.iter,prec=prec) 
        
# SCEC
thetaSCEC=EM.algorithm(y,X,beta0,Sigma0,eta0,3,distr="SCEC",max.iter=max.iter,prec=prec) 
