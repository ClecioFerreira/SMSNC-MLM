
#rm(list=ls(all=TRUE))
require(pracma)
require(numDeriv)
require(mvtnorm)
#require(sn)
require(truncdist)
require(moments)

vech <-
  function (x) {
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
      stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
  }


matrix.sqrt <- function(A)
{
   sva <- svd(A)
    if (min(sva$d)>=0)
       Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
    else
       stop("Matrix square root is not defined")
    return(Asqrt)
}

xpnd <-
  function (x, nrow = NULL) {
    dim(x) <- NULL
    if(is.null(nrow)) nrow <- (-1 + sqrt(1 + 8 * length(x))) / 2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag=TRUE)] <- 0
    output <- output + t(hold)    
    return(output)
  }



dtgamma1<-function(x,shape,rate,b=1,log=TRUE)
{
	ll<-dgamma(x,shape=shape, rate=rate,log=TRUE)-pgamma(b,shape=shape,rate=rate,log=TRUE)
	if(log==FALSE)
	{ll<-exp(ll)}
	ll
}


aux.profile.llike<-function(x,m,d,A,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	logf<-switch(dist,
			STEC=dgamma(x,shape=nu/2,rate=nu/2,log=TRUE),

			SSLEC=(log(nu)+(nu-1)*log(x))	)

	log.int<-(m/2)*log(x)-d*x/2+pcauchy(A*sqrt(x),log=TRUE)+logf
	exp(log.int)     # OK 
}

profile.llike.nu<-function(nu,y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"))
{
	p=dim(X)[1];m=dim(X)[2];n<-nrow(y)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	if(dist!="SCEC")
	{
		ll<-c()

		for(i in 1:n)
		{
			d=c(t(y[i,]-t(X[,,i])%*%beta)%*%solve(Sigma)%*%(y[i,]-t(X[,,i])%*%beta))

			A=c(t(eta)%*%(y[i,]-t(X[,,i])%*%beta))

			max.int<-switch(dist,

				STEC=Inf,

				SSLEC=1)

			aux.int<-integrate(aux.profile.llike,lower=0,upper=max.int,m=m,d=d,A=A,dist=dist,nu=nu,abs.tol=1e-15)$value

			ll[i]<-log(2)-(m/2)*log(2*pi)-(1/2)*determinant(Sigma,logarithm=TRUE)$modulus[1]+log(aux.int)
      #aux=log(2)-(m/2)*log(2*pi)-(1/2)*log(det(Sigma))+log(aux.int)
#      print(ll[i]-aux)   # OK

		}
	}
	else
	{
		nu1<-nu[1];nu2<-nu[2]

		ll<-c()

		Sigma2<-Sigma/nu2

		for(i in 1:n)
		{
			mu<-t(X[,,i])%*%beta

			A=c(t(eta)%*%(y[i,]-t(X[,,i])%*%beta))

			ll[i]<-log(2*nu1*dmvnorm(y[i,],mean=mu,sigma=Sigma2)*pcauchy(sqrt(nu2)*A)+2*(1-nu1)*dmvnorm(y[i,],mean=mu,sigma=Sigma)*pcauchy(A))

		}
	}
	-sum(ll)

}

profile.llike2<-function(theta,y,X,dist=c("STEC","SSLEC","SCEC"))
{
	p=dim(X)[1];m=dim(X)[2];n<-nrow(y)
	beta<-matrix(theta[1:p],ncol=1)
	Sigma.sqrt<-diag(c(theta[(p+1):(p+m)]))
	h<-1
	for(i in 2:m)
	{
		for(j in 1:(m-1))
		{
			if(i>j)

			{
				Sigma.sqrt[i,j]<-theta[p+m+h]

				Sigma.sqrt[j,i]<-theta[p+m+h]

				h<-h+1

			}
		}
	}
	Sigma<-Sigma.sqrt%*%Sigma.sqrt
	lambda<-matrix(theta[(p+m*(m+1)/2+1):(p+m*(m+1)/2+m)],ncol=1)
	eta=matrix.sqrt(solve(Sigma))%*%lambda
	nu<-theta[-c(1:(p+m*(m+1)/2+m))]
	if(dist!="SCEC")
	{
		ll<-c()

		for(i in 1:n)
		{
			d=c(t(y[i,]-t(X[,,i])%*%beta)%*%solve(Sigma)%*%(y[i,]-t(X[,,i])%*%beta))

			A=c(t(eta)%*%(y[i,]-t(X[,,i])%*%beta))

			max.int<-switch(dist,

				STEC=Inf,

				SSLEC=1)

			aux.int<-integrate(aux.profile.llike,lower=0,upper=max.int,m=m,d=d,A=A,dist=dist,nu=nu,abs.tol=1e-15)$value

			ll[i]<-log(2)-(m/2)*log(2*pi)-(1/2)*determinant(Sigma,logarithm=TRUE)$modulus[1]+log(aux.int)

		}
	}
	else
	{
		nu1<-nu[1];nu2<-nu[2]

		ll<-c()

		Sigma2<-Sigma/nu2

		for(i in 1:n)
		{
			mu<-t(X[,,i])%*%beta

			A=c(t(eta)%*%(y[i,]-t(X[,,i])%*%beta))

			ll[i]<-log(2*nu1*dmvnorm(y[i,],mean=mu,sigma=Sigma2)*pcauchy(sqrt(nu2)*A)+2*(1-nu1)*dmvnorm(y[i,],mean=mu,sigma=Sigma)*pcauchy(A))

		}
	}
	-sum(ll)

}

profile.llike2.sn<-function(theta,y,X,dist=c("SN","ST","SC"))
{
	p=dim(X)[1];m=dim(X)[2];n<-nrow(y)
	beta<-matrix(theta[1:p],ncol=1)
	Sigma.sqrt<-diag(c(theta[(p+1):(p+m)]))
	h<-1
	for(i in 2:m)
	{
		for(j in 1:(m-1))
		{
			if(i>j)

			{
				Sigma.sqrt[i,j]<-theta[p+m+h]

				Sigma.sqrt[j,i]<-theta[p+m+h]

				h<-h+1

			}
		}
	}
	Sigma<-Sigma.sqrt%*%Sigma.sqrt
	lambda<-matrix(theta[(p+m*(m+1)/2+1):(p+m*(m+1)/2+m)],ncol=1)
	if(dist=="ST")
	{
		nu<-theta[-c(1:(p+m*(m+1)/2+m))]

	}
	if(dist=="SN")
	{
		ll<-c()

		for(i in 1:nrow(y))
		{
			mu<-t(X[,,i])%*%beta

			ll<-dmsn(y[i,],xi=c(mu),Omega=Sigma,alpha=lambda,log=TRUE)

		}
	}
	else
	{
		if(dist=="SC")

		{nu=1}

		ll<-c()

		for(i in 1:nrow(y))
		{
			mu<-t(X[,,i])%*%beta

			ll[i]<-dmst(y[i,],xi=c(mu),Omega=Sigma,alpha=c(lambda),nu=nu,log=TRUE)

		}
	}
	-sum(ll)

}


aux.den<-function(x,y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	if(dist=="STEC" || dist=="SSLEC")
	{d=c(t(y-t(X)%*%beta)%*%solve(Sigma)%*%(y-t(X)%*%beta))}
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2,log=TRUE),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2)	)

	log.int<-pcauchy(sqrt(x)*c(t(eta)%*%(y-t(X)%*%beta)),log=TRUE)+logf
	c(exp(log.int))
}

aux.den1<-function(x,a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{  
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)	)

	log.int<-pcauchy(sqrt(x)*a)*logf
	return(log.int)
}

aux.a<-function(x,y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	if(dist=="STEC" || dist=="SSLEC")
	{d=c(t(y-t(X)%*%beta)%*%solve(Sigma)%*%(y-t(X)%*%beta))}
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2,log=TRUE),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2)	)

	log.int<-log(x)+pcauchy(sqrt(x)*c(t(eta)%*%(y-t(X)%*%beta)),log=TRUE)+logf
	c(exp(log.int))
}

aux.a1<-function(x,a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)	)

	log.int<-x*pcauchy(sqrt(x)*a)*logf
	return(log.int)
}

aux.b<-function(x,y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	if(dist=="STEC" || dist=="SSLEC")
	{d=c(t(y-t(X)%*%beta)%*%solve(Sigma)%*%(y-t(X)%*%beta))}
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2,log=TRUE),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2)	)

	log.int<-log(x)+pt(c(sqrt(3*x))*c(t(eta)%*%(y-t(X)%*%beta)),df=3,log=TRUE)+logf
	c(exp(log.int))
}

aux.b1<-function(x,a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2, log=FALSE)	)

	log.int<-x*pt(sqrt(3*x)*a,df=3)*logf
	return(log.int)
}

aux.c<-function(x,y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	if(dist=="STEC" || dist=="SSLEC")
	{d=c(t(y-t(X)%*%beta)%*%solve(Sigma)%*%(y-t(X)%*%beta))}
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2,log=TRUE),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2)	)

	log.int<-(1/2)*log(x)+dcauchy(c(sqrt(x))*c(t(eta)%*%(y-t(X)%*%beta)),log=TRUE)+logf
	c(exp(log.int))
}

aux.c1<-function(x,a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2, log=FALSE)	)

	log.int<-sqrt(x)*dcauchy(sqrt(x)*a)*logf
	return(log.int)
}

compute.a<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	if(dist=="STEC" || dist=="SSLEC")
	{
		d=c(t(y-t(X)%*%beta)%*%solve(Sigma)%*%(y-t(X)%*%beta))

		max.int<-switch(dist,

				STEC=Inf,

			SSLEC=1)

	}
	integrate(aux.a,lower=0,upper=max.int,y=y,X=X,beta=beta,Sigma=Sigma,eta=eta,dist=dist,nu=nu,abs.tol=1e-12)$value
}


compute.a1<-function(a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	{
	max.int<-switch(dist,

				STEC=Inf,

			SSLEC=1)

	}
	integrate(aux.a1,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.b<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)

	integrate(aux.b,lower=0,upper=max.int,y=y,X=X,beta=beta,Sigma=Sigma,eta=eta,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.b1<-function(a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)

	integrate(aux.b1,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.c<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)

	integrate(aux.c,lower=0,upper=max.int,y=y,X=X,beta=beta,Sigma=Sigma,eta=eta,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.c1<-function(a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)

	integrate(aux.c1,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.den<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)

	integrate(aux.den,lower=0,upper=max.int,y=y,X=X,beta=beta,Sigma=Sigma,eta=eta,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.den1<-function(a,d,p,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	max.int<-switch(dist,
			STEC=Inf,

			SSLEC=1)
  integrate(aux.den1,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.den.a<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	m<-length(y);p<-length(beta)
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	integrate(aux.den.a,lower=0,upper=Inf,y=y,X=X,beta=beta,Sigma=Sigma,eta=eta,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.abc.SCEC<-function(y,X,beta,Sigma,eta,nu)
{ 
  m=length(y) 
  p=length(beta)
	nu1<-nu[1];nu2<-nu[2]
	eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
	y<-matrix(y,ncol=1);X=matrix(X,nrow=p,ncol=m)
	Sigma2<-Sigma/nu2
	mu<-c(t(X)%*%beta)
	q<-nu1*dmvnorm(c(y),mean=mu,sigma=Sigma2)/(nu1*dmvnorm(c(y),mean=mu,sigma=Sigma2)+(1-nu1)*dmvnorm(c(y),mean=mu,sigma=Sigma))
	den<-pcauchy(sqrt(nu2)*c(t(eta)%*%(y-t(X)%*%beta)))*q+pcauchy(c(t(eta)%*%(y-t(X)%*%beta)))*(1-q)
	num.a<-nu2*pcauchy(sqrt(nu2)*c(t(eta)%*%(y-t(X)%*%beta)))*q+pcauchy(c(t(eta)%*%(y-t(X)%*%beta)))*(1-q)
	num.b<-nu2*pt(sqrt(3*nu2)*c(t(eta)%*%(y-t(X)%*%beta)),df=3)*q+pt(sqrt(3)*c(t(eta)%*%(y-t(X)%*%beta)),df=3)*(1-q)
	num.c<-sqrt(nu2)*dcauchy(sqrt(nu2)*c(t(eta)%*%(y-t(X)%*%beta)))*q+dcauchy(c(t(eta)%*%(y-t(X)%*%beta)))*(1-q)
	a<-num.a/den
	b<-num.b/den
	c0<-t(eta)%*%(y-t(X)%*%beta)*b+num.c/den
	return(list(a.theta=a,b.theta=b,c.theta=c(c0)))
}

E.step<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	n<-nrow(y) ; p=ncol(y)
	a.theta<-c();b.theta<-c();c.theta<-c()
	if(dist!="SCEC")
	{
		for(i in 1:n)
		{ # Clecio
      ai= as.numeric(t(eta)%*%(y[i,]-t(X[,,i])%*%beta))
      di= as.numeric(t(y[i,]-t(X[,,i])%*%beta)%*%solve(Sigma)%*%(y[i,]-t(X[,,i])%*%beta))  
      
      #a.a<-compute.a(y[i,],X[,,i],beta,Sigma,eta,dist=dist,nu=nu)
      a.a<-compute.a1(a=ai,d=di,p=p,dist=dist,nu=nu)

			#a.b<-compute.b(y[i,],X[,,i],beta,Sigma,eta,dist=dist,nu=nu)
      a.b<-compute.b1(a=ai,d=di,p=p,dist=dist,nu=nu)

			#a.c<-compute.c(y[i,],X[,,i],beta,Sigma,eta,dist=dist,nu=nu)
      a.c<-compute.c1(a=ai,d=di,p=p,dist=dist,nu=nu)

			#den<-compute.den(y[i,],X[,,i],beta,Sigma,eta,dist=dist,nu=nu)
      den<-compute.den1(a=ai,d=di,p=p,dist=dist,nu=nu)

			a.theta[i]<-a.a/den

			b.theta[i]<-a.b/den

			c.theta[i]<-t(eta)%*%(y[i,]-t(X[,,i])%*%beta)*b.theta[i]+a.c/den

		}
	}
	else
	{
		for(i in 1:n)
		{
			aux<-compute.abc.SCEC(y[i,],X[,,i],beta,Sigma,eta,nu)

			a.theta[i]<-aux$a.theta

			b.theta[i]<-aux$b.theta

			c.theta[i]<-aux$c.theta

		}
	}
	return(list(a.theta=a.theta,b.theta=b.theta,c.theta=c.theta))
}

M1.step<-function(y,X,Sigma,eta,a.theta,b.theta,c.theta)
{
	p=dim(X)[1];m=dim(X)[2]
	eta<-matrix(eta,ncol=1)
	sum1<-matrix(0,ncol=p,nrow=p)
	sum2<-matrix(0,ncol=p,nrow=p)
	sum3<-matrix(0,nrow=p,ncol=1)
	sum4<-matrix(0,nrow=p,ncol=1)
	sum5<-matrix(0,nrow=p,ncol=1)
	for(i in 1:nrow(y))
	{
		sum1<-sum1+a.theta[i]*X[,,i]%*%solve(Sigma)%*%t(X[,,i])

		sum2<-sum2+b.theta[i]*X[,,i]%*%(eta)%*%t(eta)%*%t(X[,,i])

		sum3<-sum3+b.theta[i]*X[,,i]%*%(eta)%*%t(eta)%*%y[i,]

		sum4<-sum4+a.theta[i]*X[,,i]%*%solve(Sigma)%*%y[i,]

		sum5<-sum5+c.theta[i]*X[,,i]%*%(eta)

	}
	solve(sum1+sum2)%*%(sum3+sum4-sum5)
}

M2.step<-function(y,X,beta,a.theta)
{
	p=dim(X)[1];m=dim(X)[2];n=nrow(y)
	beta<-matrix(beta,ncol=1)
	sum1<-matrix(0,ncol=m,nrow=m)
	for(i in 1:n)
	{
		sum1<-sum1+a.theta[i]*(y[i,]-t(X[,,i])%*%beta)%*%t(y[i,]-t(X[,,i])%*%beta)

	}
	sum1/n
}

M3.step<-function(y,X,beta,b.theta,c.theta)
{
	p=dim(X)[1];m=dim(X)[2];n=nrow(y)
	sum1<-matrix(0,ncol=m,nrow=m)
	sum2<-matrix(0,ncol=1,nrow=m)
	for(i in 1:n)
	{
		sum1<-sum1+b.theta[i]*(y[i,]-t(X[,,i])%*%beta)%*%t(y[i,]-t(X[,,i])%*%beta)

		sum2<-sum2+c.theta[i]*(y[i,]-t(X[,,i])%*%beta)

	}
	solve(sum1)%*%sum2
}

EM.algorithm<-function(y,X,beta0,Sigma0,eta0,nu0,dist=c("STEC","SSLEC","SCEC"),max.iter=1000,prec=1e-10)
{ aa=system.time({
	beta.last<-matrix(beta0,ncol=1)
	Sigma.last<-Sigma0
	eta.last<-matrix(eta0,ncol=1)
	nu.last<-nu0
  
  lower1<-switch(dist, STEC=2.001, SSLEC=1.001)
  
	i=0;dif=10
	while(i<=max.iter && dif>prec)
	{
		aux<-E.step(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)

		a.theta<-aux$a.theta

		b.theta<-aux$b.theta

		c.theta<-aux$c.theta

		beta.new<-M1.step(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)

		Sigma.new<-M2.step(y,X,beta.new,a.theta)

		eta.new<-M3.step(y,X,beta.new,b.theta,c.theta)

		if(dist!="SCEC")

		{
			nu.new<-optim(nu.last,profile.llike.nu,method="Brent",lower=lower1,upper=10,y=y,X=X,beta=beta.new,Sigma=Sigma.new,eta=eta.new,dist=dist,control = list(maxit = 20))$par
     #print(nu.new)
		}
		else

		{
			nu.new<-optim(nu.last,profile.llike.nu,method="L-BFGS-B",lower=c(0.001,0.001),upper=c(0.999,0.999),y=y,X=X,beta=beta.new,Sigma=Sigma.new,eta=eta.new,dist=dist,control = list(maxit = 40))$par

		}
		dif=abs(profile.llike.nu(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-profile.llike.nu(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))

		#dif=max(abs(c(beta.new,Sigma.new,eta.new,nu.new)-c(beta.last,Sigma.last,eta.last,nu.last)))

		eta.last=eta.new;beta.last=beta.new

		Sigma.last=Sigma.new;nu.last=nu.new

		i=i+1
   #print(i)

	}
 })
 tempo=as.numeric(aa[3])
 lognu=-profile.llike.nu(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)
 q0=ncol(Sigma.new)
 npar=length(beta.new)+length(eta.new)+q0*(q0+1)/2+length(nu.new)
 AIC=-2*lognu+2*npar
 BIC=-2*lognu+log(nrow(y))*npar
return(list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,nu=nu.new,dif=dif,iter=i,logLik=lognu,AIC=AIC,BIC=BIC,tempo=tempo))
}

EM.algorithm.fixednu<-function(y,X,beta0,Sigma0,eta0,nu,dist=c("STEC","SSLEC","SCEC"),max.iter=1000,prec=1e-10)
{
	beta.last<-matrix(beta0,ncol=1)
	Sigma.last<-Sigma0
	eta.last<-matrix(eta0,ncol=1)
	i=0;dif=10
	while(i<=max.iter && dif>prec)
	{
# print(i)
		aux<-E.step(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu)

		a.theta<-aux$a.theta

		b.theta<-aux$b.theta

		c.theta<-aux$c.theta

		beta.new<-M1.step(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)

		Sigma.new<-M2.step(y,X,beta.new,a.theta)

		eta.new<-M3.step(y,X,beta.new,b.theta,c.theta)

		dif=abs(profile.llike.nu(nu,y,X,beta.new,Sigma.new,eta.new,dist)-profile.llike.nu(nu,y,X,beta.last,Sigma.last,eta.last,dist))
     #print(profile.llike.nu(nu,y,X,beta.new,Sigma.new,eta.new,dist))
		#dif=max(abs(c(beta.new,Sigma.new,eta.new,nu.new)-c(beta.last,Sigma.last,eta.last,nu.last)))

		eta.last=eta.new;beta.last=beta.new;Sigma.last=Sigma.new
   #print(beta.new)
   #print(Sigma.new)
   #print(eta.new)

		i=i+1

	}
 lognu=profile.llike.nu(nu,y,X,beta.new,Sigma.new,eta.new,dist)
 q0=ncol(Sigma.new)
 npar=length(beta.new)+length(eta.new)+q0*(q0+1)/2+length(nu.new)
 AIC=-2*lognu+2*npar
 BIC=-2*lognu+log(nrow(y))*npar
return(list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,nu=nu,dif=dif,iter=i,logLik=lognu,AIC=AIC,BIC=BIC))
}

#EM.algorithm<-function(y,X,beta0,Sigma0,eta0,nu.min,nu.max,nu.by,dist=c("STEC","SSLEC","SCEC"),max.iter=1000,prec=1e-4)
#{
#	p=dim(X)[1];m=dim(X)[2];n=nrow(y)
#	while(nu.by>0.001)
#	{
#		nu.seq<-seq(nu.min,nu.max,nu.by)
#		res<-c()
#		for(i in 1:length(nu.seq))
#		{
#			aux<-EM.algorithm.aux(y,X,beta0,Sigma0,eta0,nu=nu.seq[i],dist=dist,max.iter=max.iter,prec=prec)
#			aux.ll<-profile.llike.nu(nu.seq[i],y,X,aux$beta,aux$Sigma,aux$eta,dist)
#			res<-rbind(res,c(aux$beta,aux$Sigma,aux$eta,aux.ll))
#		}
#		sel<-which(res[,ncol(res)]==min(res[,ncol(res)]))
#		nu.min<-nu.seq[sel]-9*nu.by/10
#		nu.max<-nu.seq[sel]+9*nu.by/10
#		nu.by<-nu.by/10
#		beta0<-matrix(res[sel,1:p],ncol=1)
#		Sigma0<-matrix(res[sel,(p+1):(p+m^2)],ncol=m,nrow=m)
#		eta0<-matrix(res[sel,(p+m^2+1):(p+m^2+m)],ncol=1)
#	}
#	return(list(beta=beta0,Sigma=Sigma0,eta=eta0,nu=nu.seq[sel]))	
#}

rsnc0<-function(n=1,Sigma,lambda)   # anterior
{
	m=length(lambda)
	u1<-abs(rnorm(n))
	u2<-abs(rnorm(n))
	u3<-rmvnorm(n,sigma=diag(m))
	lambda*u1*u2/sqrt(1+u1^2*c(t(lambda)%*%Sigma%*%lambda))+matrix.sqrt(diag(m)-u1^2*lambda%*%t(lambda)/(1+c(t(lambda)%*%Sigma%*%lambda)))%*%t(u3)
}

rsnc<-function(n=1,lambda) # Clecio
{
	m=length(lambda)
	u1<-abs(rnorm(n))
	u2<-abs(rnorm(n))
	u3<-rmvnorm(n,mean = rep(0, m), sigma=diag(m))  
	z=lambda*u1*u2/sqrt(1+u1^2*as.numeric(t(lambda)%*%lambda))+matrix.sqrt(diag(m)-u1^2*lambda%*%t(lambda)/(1+u1^2*as.numeric(t(lambda)%*%lambda)))%*%t(u3)
 return(z)
}

rsmsnc<-function(n=1,mu,Sigma,lambda,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	mu<-matrix(mu,ncol=1)
	lambda<-matrix(lambda,ncol=1)
#	u<-switch(dist,
#		STEC=rgamma(n, shape=nu/2, rate=nu/2),
#		SSLEC=rbeta(n, shape1=nu, shape2=1),
#		SCEC=ifelse(runif(n)<nu[1],nu[2],1))   # There is a problem here
  if (dist=="STEC") u= rgamma(n, shape=nu/2, rate=nu/2)
  if (dist=="SSLEC") u=rbeta(n, shape1=nu, shape2=1)
  if (dist=="SCEC") u= ifelse(runif(n)<nu[1],nu[2],1)   
#	mu+(1/u)^(1/2)*matrix.sqrt(Sigma)%*%rsnc(n=1,diag(ncol(Sigma)),lambda)
 y= mu+(1/u)^(1/2)*matrix.sqrt(Sigma)%*%rsnc(n=1,lambda)
 return(y)
}


EM.algorithm.iv<-function(y,X,nu.last,dist=c("STEC","SSLEC","SCEC"),max.iter=1000,prec=1e-4)
{ # Initial values
  
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q0=dim(X)[1]
  b0<-matrix(0,q0,q0)
  b1<-matrix(0,q0,1)
  for(i in 1:n){
    b0<-b0+X[,,i]%*%t(X[,,i])
    b1<-b1+X[,,i]%*%matrix(y[i,])
  }
  beta.last<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-as.vector(t(X[,,i])%*%beta.last)
  }
  Sigma.last<-cov(e)
  Sm.inv<-sqrtm(Sigma.last)$Binv
  #invB<-solve(B)
  lambda0<-as.matrix(skewness(e))
  eta.last=Sm.inv%*%lambda0
  theta0<-c(as.vector(beta.last),vech(Sigma.last),as.vector(eta.last),nu.last)
  #print(theta0) 
  #nu.last<-nu0
	i=1;dif=10
	while(i<=max.iter && dif>prec)
	{
		aux<-E.step(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)

		a.theta<-aux$a.theta

		b.theta<-aux$b.theta

		c.theta<-aux$c.theta

		beta.new<-M1.step(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)

		Sigma.new<-M2.step(y,X,beta.new,a.theta)

		eta.new<-M3.step(y,X,beta.new,b.theta,c.theta)

		if(dist!="SCEC")

		{
			auxnu<-optim(5,profile.llike.nu,method="Brent",lower=1e-5,upper=50,y=y,X=X,beta=beta.new,Sigma=Sigma.new,eta=eta.new,dist=dist)
    print(auxnu$value)
     nu.new=auxnu$par 
     print(nu.new)
		}
		else

		{
			nu.new<-optim(nu.last,profile.llike.nu,method="L-BFGS-B",lower=c(0.00001,0.00001),upper=c(0.99999,0.99999),y=y,X=X,beta=beta.new,Sigma=Sigma.new,eta=eta.new,dist=dist)$par

		}
		theta<-c(as.vector(beta.new),vech(Sigma.new),as.vector(eta.new),nu.new) 
    #dif=abs(profile.llike.nu(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-profile.llike.nu(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))

		#dif=max(abs(c(beta.new,Sigma.new,eta.new,nu.new)-c(beta.last,Sigma.last,eta.last,nu.last)))

    dif=sum((theta-theta0)^2)

		eta.last=eta.new;beta.last=beta.new

		Sigma.last=Sigma.new;nu.last=nu.new

		i=i+1
    theta0=theta
	}
	return(list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,nu=nu.new,dif=dif,iter=i))
}

 # Observed Fisher's Information Matrix
smsnc.MIobs<- function(y,X,theta,dist=c("STEC","SSLEC","SCEC"))
{
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,nu=nu,dif=dif,iter=i)
  p=ncol(y)
  n=nrow(y)
  beta=theta$beta
  q0=length(beta)
  #mu=X%*%beta
  #res=as.matrix(y-mu)
  Sigma=theta$Sigma
  Sinv=solve(Sigma) # Sigma^{-1}
  C=matrix.sqrt(Sigma)
  Cinv=solve(C)# Sigma^{-1/2}
  eta=theta$eta
  lambda=C%*%eta
  #Delta=Sminv%*%lambda  # eta
  q1=p*(p+1)/2
  pth=q0+p+q1
  nu=theta$nu

  dlogS=matrix(0,pth,1)
        # dlog Sigma
      for (r in 1:q1){
        indr=rep(0,q1)
        indr[r]=1
        Cr=xpnd(indr)
        auxi=Sinv%*%(Cr%*%C+C%*%Cr)
        dlogS[q0+r]=sum(diag(auxi))
        }

  dAi=matrix(0,pth,1)
  ddi=matrix(0,pth,1)
  MI=matrix(0,pth,pth)
  for (i in 1:n){
      resi=y[i,]-t(X[,,i])%*%beta  # p x 1
      ai= as.numeric(t(eta)%*%resi)
      di= as.numeric(t(resi)%*%Sinv%*%resi)
      if(dist!="SCEC")
	    {
      Iphi.i<-Iphi(a=ai,d=di,p=p,w=(p+1)/2,dist=dist,nu=nu)
      IPhi.i<-IPhi(a=ai,d=di,p=p,w=(p+2)/2,dist=dist,nu=nu)
      Ki<-IPhi(a=ai,d=di,p=p,w=p/2,dist=dist,nu=nu)
      }
      else
	    {  #print(nu)
         nu1=nu[1]
         nu2=nu[2]
         Iphi.i<-nu1*(nu2^((p+1)/2))*exp(-0.5*nu2^0.5*di)*dcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*dcauchy(ai)
         IPhi.i<-nu1*(nu2^((p+2)/2))*exp(-0.5*nu2^0.5*di)*pcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*pcauchy(ai)
         Ki<-nu1*(nu2^(p/2))*exp(-0.5*nu2^0.5*di)*pcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*pcauchy(ai)
      }
      #

      dAi[1:q0]=X[,,i]%*%Cinv%*%lambda
      dAi[(q0+q1+1):pth]=Cinv%*%resi
      ddi[1:q0]=-2*X[,,i]%*%Cinv%*%Cinv%*%resi
      for (r in 1:q1){
           indr=rep(0,q1)
           indr[r]=1
           Cr=xpnd(indr)
           ddi[q0+r]=-t(resi)%*%Cinv%*%(Cr%*%C+C%*%Cr)%*%Cinv%*%resi
           dAi[q0+r]=-t(lambda)%*%Cinv%*%Cr%*%Cinv%*%resi
           }
      scoretheta=-0.5*dlogS +(Iphi.i*dAi-0.5*IPhi.i*ddi)/Ki
      MI=MI+scoretheta%*%t(scoretheta)
    } 
 return(MI)
}

aux.Iphi<-function(x,a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)	)

	log.int<-x^w*exp(-0.5*x*d)*dcauchy(sqrt(x)*a)*logf
	return(log.int)
}

Iphi<-function(a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	{
	max.int<-switch(dist,

				STEC=Inf,

			SSLEC=1)

	}
	integrate(aux.Iphi,lower=0,upper=max.int,a,d,p,w,dist=dist,nu=nu,abs.tol=1e-12)$value
}

aux.IPhi<-function(x,a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	logf<-switch(dist,
			STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),

			SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)	)

	log.int<-x^w*exp(-0.5*x*d)*pcauchy(sqrt(x)*a)*logf
	return(log.int)
}

IPhi<-function(a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
	if(dist=="STEC" || dist=="SSLEC")
	{
	max.int<-switch(dist,

				STEC=Inf,

			SSLEC=1)

	}
	integrate(aux.IPhi,lower=0,upper=max.int,a,d,p,w,dist=dist,nu=nu,abs.tol=1e-12)$value
}


#####################################################################################################################

## SNC
EM.SNC<-function(y,X,max.iter=1000,prec=1e-5)
{ 
aa=system.time({
m=dim(X)[1]
  p=dim(X)[2]
  n=dim(X)[3]
  b0<-matrix(0,m,m)
  b1<-matrix(0,m,1)
  for(i in 1:n){
    b0<-b0+X[,,i]%*%t(X[,,i])
    b1<-b1+X[,,i]%*%matrix(y[i,])
  }
  beta.last<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-as.vector(t(X[,,i])%*%beta.last)
  }
  Sigma.last<-cov(e)
  Sm.inv<-sqrtm(Sigma.last)$Binv
  lambda.last<-as.matrix(skewness(e))
  eta.last=Sm.inv%*%lambda.last
  
  log0<- snc.loglike(y,X,beta.last,Sigma.last,lambda.last)

	iter=0;dif=1
	while(iter<=max.iter && dif>prec)
	{
   # E-step
		#aux<-E.stepSNC(y,X,c(beta.last),Sigma.last,c(eta.last))
		a.theta<-rep(1,n)                               
		b.theta<-c()                           
		c.theta<-c() 
    for (i in 1:n){
         ai= as.numeric(t(eta.last)%*%(y[i,]-t(X[,,i])%*%beta.last))
         b.theta[i]<-pt(sqrt(3)*ai,df=3)/pcauchy(ai)
         c.theta[i]<-ai*b.theta[i]+dcauchy(ai)/pcauchy(ai)
    }                        
   # M-step                           
		beta.new<-M1.step(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta) 
		Sigma.new<-M2.step(y,X,beta.new,a.theta)  
		eta.new<-M3.step(y,X,beta.new,b.theta,c.theta) 
    lambda.new<-sqrtm(Sigma.new)$B%*%eta.new
    lognu<-snc.loglike(y,X,beta.new,Sigma.new,lambda.new)
		#dif=sum(abs(c(beta.new,vech(Sigma.new),eta.new)-c(beta.last,vech(Sigma.last),eta.last))^2)
    dif<-abs(lognu-log0)
		eta.last=eta.new;beta.last=beta.new;Sigma.last=Sigma.new
		iter=iter+1
    log0=lognu 
	}
 })
 tempo=as.numeric(aa[3])
 npar=length(beta.new)+length(eta.new)+m*(m+1)/2
 AIC=-2*lognu+2*npar
 BIC=-2*lognu+log(n)*npar
return(list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,dif=dif,iter=iter,logLik=lognu,AIC=AIC,BIC=BIC,tempo=tempo))
}


# Observed Fisher's Information Matrix
snc.MIobs<- function(y,X,theta)
{
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,dif=dif,iter=i)
  p=ncol(y)
  n=nrow(y)
  beta=theta$beta
  q0=length(beta)
  Sigma=theta$Sigma
  Sinv=solve(Sigma) # Sigma^{-1}
  C=matrix.sqrt(Sigma)
  Cinv=solve(C)# Sigma^{-1/2}
  eta=theta$eta
  lambda=C%*%eta
  q1=p*(p+1)/2
  pth=q0+p+q1

  dlogS=matrix(0,pth,1)
        # dlog Sigma
      for (r in 1:q1){
        indr=rep(0,q1)
        indr[r]=1
        Cr=xpnd(indr)
        auxi=Sinv%*%(Cr%*%C+C%*%Cr)
        dlogS[q0+r]=sum(diag(auxi))
        }

  dAi=matrix(0,pth,1)
  ddi=matrix(0,pth,1)
  MI=matrix(0,pth,pth)
  for (i in 1:n){
      resi=y[i,]-t(X[,,i])%*%beta  # p x 1
      ai= as.numeric(t(eta)%*%resi)
      di= as.numeric(t(resi)%*%Sinv%*%resi)
      Iphi.i<-exp(-0.5*di)*dcauchy(ai)
      Ki<-exp(-0.5*di)*pcauchy(ai)
      
      dAi[1:q0]=X[,,i]%*%Cinv%*%lambda
      dAi[(q0+q1+1):pth]=Cinv%*%resi
      ddi[1:q0]=-2*X[,,i]%*%Cinv%*%Cinv%*%resi
      for (r in 1:q1){
           indr=rep(0,q1)
           indr[r]=1
           Cr=xpnd(indr)
           ddi[q0+r]=-t(resi)%*%Cinv%*%(Cr%*%C+C%*%Cr)%*%Cinv%*%resi
           dAi[q0+r]=-t(lambda)%*%Cinv%*%Cr%*%Cinv%*%resi
           }
      scoretheta=-0.5*dlogS +(Iphi.i*dAi-0.5*Ki*ddi)/Ki
      MI=MI+scoretheta%*%t(scoretheta)
    } 
 return(MI)
}

snc.loglike<- function(y,X,beta0,Sigma,eta)
{
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,dif=dif,iter=i)
  n=nrow(y)
  logv=0
  for (i in 1:n){
      mui<-t(X[,,i])%*%beta0
      ai= as.numeric(t(eta)%*%(y[i,]-mui))
      aux<-2*dmvnorm(y[i,],mean=as.vector(mui),sigma=Sigma)*pcauchy(ai)
      logv=logv+log(aux)
    } 
 return(logv)
}

snc.log<- function(y,X,theta)
{
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,dif=dif,iter=i)
  n=nrow(y)
  beta0=theta$beta
  Sigma=theta$Sigma
  eta=theta$eta
  logv=0
  for (i in 1:n){
      mui<-t(X[,,i])%*%beta0
      ai= as.numeric(t(eta)%*%(y[i,]-mui))
      aux<-2*dmvnorm(y[i,],mean=as.vector(mui),sigma=Sigma)*pcauchy(ai)
      logv=logv+log(aux)
    } 
 return(logv)
}

