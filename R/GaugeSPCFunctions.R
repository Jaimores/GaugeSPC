#######                              Gauge - SPC                                  #######
#                             Jaime Mosquera Restrepo                             #######
#######                              2016 -2019                                   #######


#####Summary  #######################################################################################
#This scripts contains the  function's set employed for make a performance's evaluation and design  # 
#of proposed based gauge Control Charts                                                             #
#****************************************************************************************************

library(sn)   #Required for skew-normal distibutions
library(GA)   #Required for the Genetic Algortihm
library(SparseM) #Required for build the matrix image in Markov Chain
library(svMisc)  #Visualitation of progress fucntion
#######Preliminaries. Required Functions ############################################################
# Set of auxiliar function required by the principal function or for comparing with alternative CC  #
#****************************************************************************************************

#Find the shape parameter for Lognormal, Skewnormal and Weibull distribution with a prefixed skewness  
Root.skew<-function(skew,dist,tol=0.000001,maxiter=100){ 

  k=0
  
  if(dist=="lognorm"){
  if(skew<=0){stop("for Lognormal distribution, skewness must be positive")}
  
  
  f<-function(x,skew){
    f1<-2*x+log(exp(x)+3)-log(skew^2+4)
    fd<-2+exp(x)/(exp(x)+3)
    return(c(f1,fd))
  }
  
  xk=0.5
  repeat{
    p=f(xk,skew)
    dx=p[1]/p[2]
    xk1=xk-dx
    xk=xk1
    k=k+1
    if (abs(dx)<=tol|k>=maxiter){break}
  }
  
  }
  
  else if(dist=="weibull"){
  
    if(skew<(-1.1)){stop("for Weibull distribution, skewness must be upper than -1.1")}  
  
  f<-function(x,skew){
    gam1=gamma(1+(1/x));gam2=gamma(1+(2/x));gam3=gamma(1+(3/x)) 
    phi1=digamma(1+(1/x));phi2=digamma(1+(1/x));phi3=digamma(1+(1/x))
    
    f1<-(gam3-3*gam1*gam2+2*(gam1^3)-skew*((gam2-(gam1^2))^(3/2)))
    fd<-(-3/(x^2))*(gam3*phi3-gam1*gam2*(phi1-2*phi2)+(gam1^3)*phi1-skew*(gam2*phi2-(gam1^2)*phi1)*sqrt(gam2-(gam1^2)))
    return(c(f1,fd))
    }  
  
  if((skew >=0.5)){a<-0.1;b<-2.3}
  else if((skew<0.5)&(skew>=-0.5)){a<-2.0;b<-10}
  else if((skew<(-0.5))&(skew>=(-1))){a<-8;b<-50}
  else if((skew<(-1))&(skew>=(-1.1))){a<-40;b<-160}
  
  if(skew<(-1)){tol=0.0000000001}
    repeat{
      c<-(a+b)/2
      fc<-f(c,skew)[1]
      if(fc*f(b,skew)[1]>0){b<-c}else{a<-c}
      k<-k+1
      if (abs(fc)<tol|k>maxiter){break}
    } 
  xk=c
  }
  
  if(dist=="lognorm"){xsal=sqrt(xk);names(xsal)="sigmalog"}
  else if(dist=="weibull"){xsal=xk;names(xsal)="a"}
  
  return(c(xsal))
}


#Lognormal 3-parameters distribution (gamma, meanlog, sdlog)

dlnorm3<-function(x,gamma,meanlog,sdlog){
  x<-(x-gamma)
  den<-dlnorm(x,meanlog,sdlog)
  return(den)
}

plnorm3<-function(x,gamma,meanlog,sdlog){
  x<-(x-gamma)
  prob<-plnorm(x,meanlog,sdlog)
  return(prob)
}

qlnorm3<-function(q,gamma,meanlog,sdlog){
  xp<-gamma + qlnorm(q,meanlog,sdlog)
  return(xp)
}

rlnorm3<-function(n,gamma,meanlog,sdlog){
  s<-gamma+rlnorm(n,meanlog,sdlog)
  return(s)
}

#Tails probability and parameters associated to (Normal, Lognormal2, Lognormal3, Weibull3, Skewnormal) distributions
# when gauges dimensions are allocated in (q) ang there are a shift (delta,r) in the proccess. 

conv<-function(delta,r,q,mu=10,sig=1,dist="norm",skew=0){ #mu y sig en escala original
  
  if (q<=0|q>=1) {stop("Parametros Inconsistentes")}
  
  if(dist=="norm"){
    z=-qnorm(q/2)
    qsd= pnorm((-z-delta)/r)
    qld= 1-pnorm((z-delta)/r)
    S=mu-z*sig
    L=mu+z*sig
  }
  
  else if(dist=="weibull"){
    
    a=Root.skew(skew,"weibull")
    gam1=gamma(1+(1/a));gam2=gamma(1+(2/a));
    b=sig/sqrt(gam2-(gam1^2))
    gamma=mu-b*gam1
    S=gamma+b*((-log(1-q/2))^(1/a))
    L=gamma+b*((-log(q/2))^(1/a))
    
    
    b1=r*b
    gamma1=gamma+b*(gam1*(1-r)+delta*sqrt(gam2-(gam1^2)))
    
    fact=gam1*(1-r)+delta*sqrt(gam2-(gam1^2))

    qsd=pweibull(S-gamma1,a,b1)
    qld=1-pweibull(L-gamma1,a,b1)
    #qsd= 1-exp((-1/(r^a))*((-log(1-q/2))^(1/a)-fact)^a)
    #qld=exp((-1/(r^a))*((-log(q/2))^(1/a)-fact)^a)
  }
  
  
  else if(dist=="lognorm2"){
    mulog=log(mu/(sqrt(1+(sig^2)/(mu^2))))
    siglog=sqrt(log(1+ (sig^2)/(mu^2)))
    S=qlnorm(q/2, mulog, siglog)
    L=qlnorm(1-q/2, mulog, siglog)
    mu1=mu+delta*sig; sig1=r*sig
    mulog1=log(mu1/(sqrt(1+(sig1^2)/(mu1^2))))
    siglog1=sqrt(log(1+ (sig1^2)/(mu1^2)))
    
    qsd= plnorm(S, mulog1, siglog1)
    qld= 1-plnorm(L, mulog1, siglog1)  
    gamma=gamma1=0
  }
  else if (dist=="lognorm3"){
    z=-qnorm(q/2)
    siglog=Root.skew(skew,"lognorm")
    mulog=log(sig/sqrt(exp(siglog^2)-1))-0.5*(siglog^2)
    gamma=mu-sig/sqrt(exp(siglog^2)-1)
    
    siglog1=siglog
    mulog1=log(r)+mulog
    gamma1=gamma +exp(mulog+0.5*(siglog^2))*(1-r+delta*sqrt(exp(siglog^2)-1))
    
    S=gamma+exp(mulog-siglog*z)
    L=gamma+exp(mulog+siglog*z)
    #qsd= plnorm3(S,gamma1, mulog1, siglog1)
    #qld= 1-plnorm3(L,gamma1, mulog1, siglog1)
    
    
    qsd=pnorm(-z+(1/siglog)*log(max(0,(1/r)+exp(0.5*(siglog^2)+z*siglog)*((1-1/r)-delta*sqrt(exp(siglog^2)-1)/r))))
    qld=1-pnorm(z+(1/siglog)*log(max(0,(1/r)+exp(0.5*(siglog^2)-z*siglog)*((1-1/r)-delta*sqrt(exp(siglog^2)-1)/r)))) # maximo para evitar que el log se indetermine con valores extremos de q
  }
  else if (dist=="snorm"){
    d=sign(skew)*sqrt((pi/2)*(abs(skew)^(2/3))/(abs(skew)^(2/3) + ((4-pi)/2)^(2/3)))
    lambda=d/sqrt(1-d^2)
    zs=qsn(q/2,0,1,lambda)
    zl=qsn(1-q/2,0,1,lambda)
    omega=sig/sqrt(1-2*(d^2)/pi)
    xi=mu-omega*d*sqrt(2/pi)
    S=xi+zs*omega
    L=xi+zl*omega
      
    qsd= psn((zs-sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2)))/r,0,1,lambda)
    qld= 1-psn((zl-sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2)))/r,0,1,lambda)
    
    omega1=r*omega
    xi1=xi+omega*sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2))
  }
  
  else{stop("dist must be: norm, snorm, lognorm2, lognorm3, weibull")}
  
  if(dist=="norm"){sol<-c(qsd,qld,S,L)
  names(sol)<-c("qsd","qld","S","L")}
  else if(dist=="snorm"){sol<-c(qsd,qld,S,L,xi,omega,lambda,xi1,omega1,lambda)
  names(sol)<-c("qsd","qld","S","L","xi","omega","lambda","xi1","omega1","lambda1")
  }
  else if(dist=="weibull"){sol<-c(qsd,qld,S,L,gamma,b,a,gamma1,b1,a)
  names(sol)<-c("qsd","qld","S","L","gamma","b","a","gamma1","b1","a1")
  }
  else{ sol<-c(qsd,qld,S,L,gamma,mulog,siglog,gamma1,mulog1,siglog1)
   names(sol)<-c("qsd","qld","S","L","gamma","mulog","siglog","gamma1","mulog1","siglog1")
  }
   return(sol)
}

##Mean, deviation and skewnees for a (weibull, Skewnormal and Lognormal3) distribution with parameters (location, escale, shape).
ref.weibull<-function(gamma,b,a){
  gam1=gamma(1+(1/a));gam2=gamma(1+(2/a));gam3=gamma(1+(3/a))
  mean<-gamma+b*gam1
  sd<-b*sqrt(gam2-(gam1^2))
  CA<-round(((gam3-3*gam1*gam2+2*(gam1^3))/((gam2-(gam1^2))^(3/2))) ,3)     
  sol<-c(mean,sd,CA)
  names(sol)<-c("mean","sd","CA")
  return(sol)
}

ref.snorm<-function(xi,omega,lambda){
  d=lambda/sqrt(1+lambda^2)
  mean<-round(xi+omega*d*sqrt(2/pi),2)
  sd<-round(omega*sqrt(1-2*(d^2)/pi),2)
  CA<-round(((4-pi)/2)*(d*sqrt(2/pi)/sqrt(1-2*(d^2)/pi))^3,3)     
  sol<-c(mean,sd,CA)
  names(sol)<-c("mean","sd","CA")
  return(sol)
}

ref.lognorm3<-function(gamma,mulog,siglog){
  mean<-round(gamma+exp(mulog+0.5*(siglog^2)),2)
  sd<-round(sqrt(exp(2*mulog+(siglog^2))*(exp(siglog^2)-1)),2)
  CA<-round((exp(siglog^2)+2)*sqrt(abs(exp(siglog^2)-1)),3)
  sol<-c(mean,sd,CA)
  names(sol)<-c("mean","sd","CA")
  return(sol)
}


#Cumulative Probability in Omega-Incontrol for wYSYL

Fmultw<-function(n,Lim,w,qs,ql){
Prob=0
if ((min(qs,ql) <= 0.0)|((qs + ql) > 1.0)| (n <= 0)| (Lim<=0)) {stop("Parametros Inconsistentes")}
else {
   for (s in 0:n){
      for (l in 0:(n-s)){
         if (((w*s+l)<Lim)&((s+w*l)<Lim)){Prob = Prob + (factorial(n)/(factorial(s)*factorial(l)*factorial(n-s-l)))*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))}
                         }
                 }
}
return(min(Prob,0.999999999999999999999999999999))
}

# Probability distribution function of wYsYl = max(wYs + Yl; Ys + wYl)

pwYsYl<-function(w,n,qs,ql){
  
  if ((min(qs,ql) < 0.0)|((qs + ql) > 1.0)| (n <= 0)) {stop("Parametros Inconsistentes")}
  else {
    wYSYL=Prob=0
    i=1
    for (s in 0:n){
      for (l in 0:(n-s)){
        wYSYL[i]<-max(w*s+l,s+w*l)
        Prob[i] = (factorial(n)/(factorial(s)*factorial(l)*factorial(n-s-l)))*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))
        i=i+1
      }
    }
  }
  est<-sort(unique(wYSYL))
  p<-0;Facum<-0
  acum<-0
  
  
  for(j in 1:length(est)){
    p[j]<-sum(Prob[wYSYL==est[j]])
    Facum[j]<-acum+p[j]
    acum<-Facum[j] 
  }
  salida<-cbind(est,p,Facum)
  colnames(salida)<-c("wYSYL","Prob","Facum")
  return(salida)  
}




###Equivalences with Steiner Approach#########################################

# What are the equivalents (q,w) of Steiner's optimal proposal?
wqSte2<-function(delta,r,opt){
   if(opt=="mean"){q=2*pnorm(-0.612)}else if (opt=="sd"){q=2*pnorm(-1.4825)}else if(opt=="mix"){q=2*pnorm(-0.8487)}else{stop("Invalid Scheme Opt= (mean,sd,mix)")}
   qd<-conv(delta,r,q)[1:2]
   p1=log(q/2/qd[1]);p2=log((1-q)/(1-sum(qd)));p3=log(q/2/qd[2])
   w=(p1-p2)/(p3-p2)
   sol<-c(w,q)
   names(sol)<-c("w","q")
   return(sol)
}

# Fixed (q), What is the equivalent (w) of Steiner's proposal?
wSte2<-function(delta,r,q){
   qd<-conv(delta,r,q)[1:2]
   p1=log(q/2/qd[1]);p2=log((1-q)/(1-sum(qd)));p3=log(q/2/qd[2])
   w=(p1-p2)/(p3-p2)
   names(w)<-"w"
   return(w)
}

### Performance of Shewhart Control Charts Xb, S, Xb,S###########################


#ARL for Shewhart Control Chart
ARLXb<-function(n,alp,delta,r){

  Z=qnorm(alp/2,lower.tail = F)  
  beta=pnorm((Z-delta*sqrt(n))/r)-pnorm((-Z-delta*sqrt(n))/r) 
  ARL=1/(1-beta)
  return(ARL)   
}

ARLS<-function(n,alp,r){
  chi<-qchisq(alp,n-1,lower.tail=F)
  ARL<-1/pchisq(chi/(r^2),n-1,lower.tail=F)
  return(ARL)   
}

ARLXS<-function(n,alp,delta,r){
  
  alpc<-1-sqrt(1-alp); Z<-qnorm(alpc/2,lower.tail=F);chi<-qchisq(alpc,n-1,lower.tail=F)   
  
  betaX<-pnorm((Z-delta*sqrt(n))/r)-pnorm((-Z-delta*sqrt(n))/r)
  betaS<-pchisq(chi/(r^2),n-1,lower.tail=T)
  ARL<-1/(1-betaX*betaS)
  return(ARL)   
}

ARLXb.TMV<-function(n1,n2,L,w,delta,r){
  ARL=0
  ARL[1]= 1/(2*pnorm(-L))
  
  p1= pnorm(w)-pnorm(-w)
  p2= 1-2*pnorm(-L)-p1
  
  p1= p1/(p1+p2)
  p2= 1-p1
  ARL[2]=  p1*n1+(1-p1)*n2 
  
  p11= pnorm((w-delta*sqrt(n1))/r)-pnorm((-w-delta*sqrt(n1))/r)
  p12= pnorm((L-delta*sqrt(n1))/r)-pnorm((-L-delta*sqrt(n1))/r)-p11
  p21= pnorm((w-delta*sqrt(n2))/r)-pnorm((-w-delta*sqrt(n2))/r)
  p22= pnorm((L-delta*sqrt(n2))/r)-pnorm((-L-delta*sqrt(n2))/r)-p21
  
  ARL[3]= (p1*(1-p22+p12)+p2*(1-p11+p21))/((1-p11)*(1-p22)-p21*p12) 
  
  return(round(ARL,3))
}

#Curve ARL for Shewhart Control Chart

CARLXb<-function(n,alp,delta,r,curve="mean"){
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
      ARL[i]<-ARLXb(n,alp,0,r[i])
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
      qd<-conv(delta[i],1,q,mu,sig,gamma,dist)[1:2]
      ARL[i]<-ARLXb(n,alp,delta[i],1)
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
}

CARLS<-function(n,alp,r){
  ARL<-ARLS(n,alp,r)  
  return(cbind(r,ARL))
}

CARLXS<-function(n,alp,delta,r,curve="mean"){
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
      ARL[i]<-ARLXS(n,alp,0,r[i])
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
      qd<-conv(delta[i],1,q,mu,sig,gamma,dist)[1:2]
      ARL[i]<-ARLXS(n,alp,delta[i],1)
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
  
}

CARLXb.TMV<-function(n1,n2,L,w,delta,r,curve="mean"){
  
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
      ARL[i]<-ARLXb.TMV(n1,n2,L,w,0,r[i])
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
      qd<-conv(delta[i],1,q,mu,sig,gamma,dist)[1:2]
      ARL[i]<-ARLXb.TMV(n1,n2,L,w,delta[i],1)
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
}


#####Block 1 - Control Chart for Univariate Gauge - Fixed Sample#####################################
##Functions used to evaluate the performance of proposals for univariate control gauge              # 
#****************************************************************************************************


#ARl evaluation for univariate gauge control chart 

ARLwYsYl<-function(n,Lc,w,qs,ql){
   Prob=0
   if ((min(qs,ql) < 0.0)|((qs + ql) > 1.0)| (n <= 0)| (Lc<=0)) {(ARL=NA)}
   else {
     fact.n<-factorial(n)
      for (s in 0:n){
        fact.s<-factorial(s)
         for (l in 0:(n-s)){
               if (((w*s+l)<Lc)&((s+w*l)<Lc)){Prob = Prob + (fact.n/(fact.s*factorial(l)*factorial(n-s-l)))*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))}
         }
      }
     ARL=1/(1-round(Prob,8))
   }
   return(ARL) }

#Find the Steiner approach parameters equivalents a wYSYL and nearly to ARLdes
parSte2<-function(n,delta,r,opt,ARLdes){
  if((delta==0)&(r<=1)){stop("delta and r must be different to (0,1)")}
  
  par<-wqSte2(delta,r,opt)
  k=1;i=0;lim=0;ARL=0
  
  for (s in 0:n){
    for (l in 0:(n-s)){
      lim[k]<-s+par[1]*l;k=k+1
    }
  }
  lim<-sort(lim,decreasing=T)
  
  repeat{
    i=i+1
    if (ARLwYsYl(n,lim[i],par[1],par[2]/2,par[2]/2)<ARLdes|i==length(lim)) break
  }
  
  if(i>1){ind<-c(i-1,i)}else{ind<-1}
  sol<-c(lim[ind],par)
  names(sol)<-c("CL+","CL-","w","q")
  return(sol)
}

#ARl for joint YS-YL, YT Control Chart - any probability distibution
ARLDSYsYl<-function(n,LcD,LcS,qs,ql){
  Prob=0
  if ((min(qs,ql) < 0.0)|((qs + ql) > 1.0)| (n <= 0)| (LcD<=0)| (LcS<=0)) {(ARL=NA)}
  else {
    fact.n<-factorial(n)
    for (s in 0:n){
      fact.s<-factorial(s)
      for (l in 0:(n-s)){
        if (((-s+l)<LcD)&((s-l)<LcD)&((s+l)<LcS)){Prob = Prob + (fact.n/(fact.s*factorial(l)*factorial(n-s-l)))*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))}
      }
    }
    ARL=1/(1-round(Prob,8))
  }
  return(ARL)                         }

#Distance of Ys,YL from  E[YS],E[YL]-Required for D2YsYl CC

D2YsYl<-function(Ys,Yl,n,q,origin="expected"){
  if(origin=="zero"){ mu<-0}else{ mu<-n*(q/2)}
  ro<-(-(q/2)/(1-q/2));sig2<-n*(q/2)*(1-q/2)
  D2<-(1/((1-(ro^2))*(sig2)))*((Ys-mu)^2+(Yl-mu)^2-2*ro*(Ys-mu)*(Yl-mu))
  return(D2)                 }

#ARL for D2YSYL Control Chart - Ellipsoidal control for any probability distribution.
ARLD2YsYl<-function(n,Lc,q,qs,ql){
  Prob=0
  if ((min(qs,ql) < 0.0)|((qs + ql) > 1.0)| (n <= 0)| (Lc<=0)) {(ARL=NA)}
  else {
    fact.n<-factorial(n)
    for (s in 0:n){
      fact.s<-factorial(s)
      for (l in 0:(n-s)){
        if (D2YsYl(s,l,n,q)<Lc){Prob = Prob + (fact.n/(fact.s*factorial(l)*factorial(n-s-l)))*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))}
      }
    }
    ARL=1/(1-min(round(Prob,8),0.99999999999))
  }
  
  return(ARL) 
  }

#Curve ARL for fixed sample gauge control chart

CARLwYsYl<-function(n,Lc,w,q,delta,r,curve="mean",mu=10,sig=1,gamma=0,dist="norm"){ 
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
      qd<-conv(0,r[i],q,mu,sig,gamma,dist)[1:2]
      ARL[i]<-ARLwYsYl(n,Lc,w,qd[1],qd[2])
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
      qd<-conv(delta[i],1,q,mu,sig,gamma,dist)[1:2]
      ARL[i]<-ARLwYsYl(n,Lc,w,qd[1],qd[2])
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
  
}

#ARL wYsYl-EWMA con numero fijo de etados (m)

ARLwYsYl.EWMA<-function(n,Lim,w,lambda,q,delta,r,m,opt=1,mu=10,sig=1,dist="norm",skew=0,image=F){
  
  #generar tabla con distribucion de wYSYL
  if((lambda > 1.0)| (lambda<=0)| (Lim<=0)) {stop("Parametros Inconsistentes")}
  
  Mst=Mlt=Mt=Prob.IC=Prob.OOC=0;p.IC=p.OOC=0;     #se inician los vectores para construir la distribución
  h=0                                             #Posicion inicial para distribucion  bajo control
  Mu<-0;E2<-0                                     #Valor Inicial para la media y varianza
  
  qd=conv(delta,r,q,mu,sig,dist,skew)[1:2]
  qs<-qd[1]; ql<-qd[2]                            # obtiene probabilidades en cola
  
  
  
  fact.n<-factorial(n)                            #para evitar repetir el calculo
  
  
  
  if(opt==1){
    # construye la distribuci?n de max(wYs+Yl, Ys+wYl)
    for (s in 0:n){
      fact.s<-factorial(s)
      for (l in 0:(n-s)){
        h=h+1
        cte<-fact.n/(fact.s*factorial(l)*factorial(n-s-l))
        Mt[h]<-max(w*s+l,s+w*l)
        Prob.IC[h] = cte*(q/2)^(s+l)*((1-q)^(n-s-l))
        Prob.OOC[h] = cte*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))
      }
    }
    
    est<-sort(unique(Mt))  
    
    #resume la distribución, eliminando valores repetidos
    for(k in 1:length(est)){
      p.IC[k]<-sum(Prob.IC[Mt==est[k]])
      Mu<-Mu+est[k]*p.IC[k]
      E2<-E2+(est[k]^2)*p.IC[k]
      p.OOC[k]<-sum(Prob.OOC[Mt==est[k]])
    }
    
    #Calcula estad?sticos resumen, ------ojo Min puede variar
    Min<-est[1]
    Var<-E2-(Mu^2)
    LC<-Mu+Lim*sqrt(Var*lambda/(2-lambda))            #dado que L>0, se garantiza que LC siempre es mayor que mu
    #Min<-Mu-5*sqrt(Var*lambda/(2-lambda)) 
    
    #Crea y llena Matrice de probabilidades de transici?n bajo control P, fuera de control Q.
    
    amp=(LC-Min)/m                                  #Amplitud de los intervalos de Zt
    P=Q=matrix(0,ncol=m,nrow=m)                       #Inicio de probabilidades de transicion bajo y fuera de control.
    
    for(i in 1:m){                                    #Recorre todos los estados iniciales
      PMi=Min+(i-0.5)*amp                             #Punto medio del estado de partida
      
      for(h1 in 1:length(est)){                       #Para todos los valores del estadistico (est), evalua el intervalo de llegada
        z=lambda*est[h1]+(1-lambda)*PMi              #Solo modifica los P[i,j] requeridos
        
        if(z<LC){
          j=1+floor(max((z-Min)/amp,0))               #Identifica el estado j al que llegar? Zt
          P[i,j]=P[i,j]-p.IC[h1]                      #Acumula probabilidad en P[i,j]
          Q[i,j]=Q[i,j]-p.OOC[h1]}                    #Acumula probabilidad en Q[i,j]
      }
    }
    i.prom=1+floor(max((Mu-Min)/amp,0))               # intervalo asociado al valor medio
    
  }else{
    # construye la distribuci?n de wYs+Yl y de Ys+wYl
    for (s in 0:n){                                 
      fact.s<-factorial(s)
      for (l in 0:(n-s)){
        cte<-fact.n/(fact.s*factorial(l)*factorial(n-s-l))
        h=h+1
        Mst[h]<-w*s+l; Mlt[h]<-s+w*l
        Prob.IC[h] =  cte*(q/2)^(s+l)*((1-q)^(n-s-l))
        Prob.OOC[h] = cte*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))
      }
    }
    
    Min<-min(0,w*n)
    Mu=n*(q/2)*(w+1)   #valor medio de ambos, Mst, Mlt
    Var=n*q/2*(1-q/2)*(w*w+1-w*q/(1-q/2))
    LC<-Mu+Lim*sqrt(Var*lambda/(2-lambda))
    
    
    #####se construyen las matrices de probabilidades de transicion
    
    mT=m*m                                         #m estados por cada estadistico
    amp=(LC-Min)/m                                #Amplitud de los intervalos de Zt
    P=Q=matrix(0,ncol=mT,nrow=mT)
    
    for(i.s in 1:m){                              # recorre todos los estados i para Ms
      PMi.s=Min+(i.s-0.5)*amp                     # Punto medio del estado i.s
      
      for(i.l in 1:m){                            # recorre todos los estados i para Ml
        PMi.l=Min+(i.l-0.5)*amp                   # Punto medio del estado i.l
        k=(i.l-1)*m+i.s                           ## transforma estados i.s,i.l a estado bidimensional k.
        
        for(h1 in 1:h){                           #Recorre todos el esapcio muestral, identificando el estado de llegada
          Zs=lambda*Mst[h1]+(1-lambda)*PMi.s
          Zl=lambda*Mlt[h1]+(1-lambda)*PMi.l
          
          if(max(Zs,Zl)<LC){
            j.s=1+floor(max((Zs-Min)/amp,0))      #Obtiene el estado de destino j.s para cada Zts
            j.l=1+floor(max((Zl-Min)/amp,0))      #Obtiene el estado de destino j.l para cada Ztl
            l=(j.l-1)*m+j.s                       #transforma j.s, j.l a Estado bidimensional de llegada
            P[k,l]=P[k,l]-Prob.IC[h1]             # Acumula las probabilidades en P[k,l]
            Q[k,l]=Q[k,l]-Prob.OOC[h1]}           # Acumula las probabilidades en Q[k,l]
        }   #acumula en Negativo la probabilidades
      } 
    }
    
    i.lprom=i.sprom=1+floor(max((Mu-Min)/amp,0))  #Intervalos asociados a valor medio de Ztl Zts 
    i.prom=(i.lprom-1)*m+i.sprom                  #Intervalo asociado a valor medio (estado bidimensional)
  }
  
  
  
  #Calculo del ARL - Cadena de Markov
  #if(opt==1){Id=diag(1,nrow=m)}else{Id=diag(1,nrow=mT)}
  if(opt==1){mT=m}else{mT=mT}
  
  for (i in 1:mT){
       P[i,i]=1+P[i,i]
       Q[i,i]=1+Q[i,i]    #Ahora P y Q contienes a Id-P e Id - Q respectivamente
  }
  if(image==T){
            x11()
            par(mfrow=c(1,2))
            image(as.matrix.csr(P))}
  if(rcond(P)<2.220446e-16){
    ARL0.ZS=10000; ARL1.ZS=10000;ARL0.SS=10000; ARL1.SS=10000 
  }else{
  P=solve(P); Q=solve(Q)}    # se reemplazan por las inversas, para evitar consumo de memoria
  
  
  if(image==T){image(as.matrix.csr(P))}
  
  ARL0.ZS=sum(P[i.prom,])
  ARL1.ZS=sum(Q[i.prom,])
  P.0=P[i.prom,]/ARL0.ZS
  
  ARL0.SS= sum(P.0%*%P)
  ARL1.SS= sum(P.0%*%Q)
  salida<-c(ARL0.ZS,ARL1.ZS,ARL0.SS,ARL1.SS)
  return(salida)
}

#division de estados segun amplitud deseada para cada estado(amp)
ARLwYsYl.EWMA<-function(n,Lim,w,lambda,q,delta,r,opt=1,mu=10,sig=1,dist="norm",skew=0,amp=0.01,image=F){
          
          #generar tabla con distribucion de wYSYL
          if((lambda > 1.0)| (lambda<=0)| (Lim<=0)) {stop("Parametros Inconsistentes")}
          
          Mst=Mlt=Mt=Prob.IC=Prob.OOC=0;p.IC=p.OOC=0;     #se inician los vectores para construir la distribución
          h=0                                             #Posicion inicial para distribucion  bajo control
          Mu<-0;E2<-0                                     #Valor Inicial para la media y varianza
          
          
          qd=conv(delta,r,q,mu,sig,dist,skew)[1:2]
          qs<-qd[1]; ql<-qd[2]                            # obtiene probabilidades en cola
          fact.n<-factorial(n)                            #para evitar repetir el calculo
          
          if(opt==1){
                    # construye la distribuci?n de max(wYs+Yl, Ys+wYl)
                    for (s in 0:n){
                              fact.s<-factorial(s)
                              for (l in 0:(n-s)){
                                        h=h+1
                                        cte<-fact.n/(fact.s*factorial(l)*factorial(n-s-l))
                                        Mt[h]<-max(w*s+l,s+w*l)
                                        Prob.IC[h] = cte*(q/2)^(s+l)*((1-q)^(n-s-l))
                                        Prob.OOC[h] = cte*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))
                              }
                    }
                    
                    est<-sort(unique(Mt))  
                    
                    #resume la distribuci?n, eliminando valores repetidos
                    for(k in 1:length(est)){
                              p.IC[k]<-sum(Prob.IC[Mt==est[k]])
                              Mu<-Mu+est[k]*p.IC[k]
                              E2<-E2+(est[k]^2)*p.IC[k]
                              p.OOC[k]<-sum(Prob.OOC[Mt==est[k]])
                    }
                    
                    #Calcula estad?sticos resumen, ------ojo Min puede variar
                    Min<-round(lambda*(est[1]+((1-lambda)^9)*(est[1]+Mu/lambda)),2) # ser? el primer punto medio  
                    Var<-E2-(Mu^2)
                    LC<-round(Mu+Lim*sqrt(Var*lambda/(2-lambda)),2)            #dado que L>0, se garantiza que LC siempre es mayor que mu
                    #Crea y llena Matrice de probabilidades de transici?n bajo control P, fuera de control Q.
                    
                    if(LC<=Min|LC<=Mu){
                              ARL0.ZS=1; ARL1.ZS=1;ARL0.SS=1; ARL1.SS=1 
                              salida<-c(ARL0.ZS,ARL1.ZS,ARL0.SS,ARL1.SS)
                              return(salida)}
                    
                    #Amplitud de los intervalos de Zt
                    m=ceiling((LC-Min)/amp)                  ###A pesar de redondear, hay problemas con la division, no es entera
                    
                    
                    P=Q=matrix(0,ncol=m,nrow=m)                       #Inicio de probabilidades de transicion bajo y fuera de control.
                    i.prom=floor(1+max((Mu-Min)/amp,0))             # intervalo asociado al valor medio
                    #
                    
                    for(i in 1:m){                                    #Recorre todos los estados iniciales
                              PMi=Min+(i-0.5)*amp                               #Punto medio del estado de partida
                              
                              for(h1 in 1:length(est)){                       #Para todos los valores del estadistico (est), evalua el intervalo de llegada
                                        z=lambda*est[h1]+(1-lambda)*PMi    #Solo modifica los P[i,j] requeridos
                                        
                                        if(z-LC<(-5e-9)){
                                                  j=floor(1+max((z-Min)/amp,0))             #Identifica el estado j al que llegar? Zt
                                                  if(j>m){
                                                            print(paste(LC,",",w,",",q,",",Lim,",",lambda))
                                                            print(paste(z,",",Min,",",Mu,",",m,",",j,",",i.prom,",",nrow(P),ncol(P)))}
                                                  P[i,j]=P[i,j]-p.IC[h1]                      #Acumula probabilidad en P[i,j]
                                                  Q[i,j]=Q[i,j]-p.OOC[h1]}                    #Acumula probabilidad en Q[i,j]
                              }
                    }
                    
                    
          }else{
                    # construye la distribuci?n de wYs+Yl y de Ys+wYl
                    for (s in 0:n){                                 
                              fact.s<-factorial(s)
                              for (l in 0:(n-s)){
                                        cte<-fact.n/(fact.s*factorial(l)*factorial(n-s-l))
                                        h=h+1
                                        Mst[h]<-w*s+l; Mlt[h]<-s+w*l
                                        Prob.IC[h] =  cte*(q/2)^(s+l)*((1-q)^(n-s-l))
                                        Prob.OOC[h] = cte*(qs^s)*(ql^l)*((1-qs-ql)^(n-s-l))
                              }
                    }
                    
                    Min<-round(lambda*(min(0,w*n)+((1-lambda)^1)*(min(0,w*n)+Mu/lambda)),2)
                    Mu=n*(q/2)*(w+1)   #valor medio de ambos, Mst, Mlt
                    Var=n*q/2*(1-q/2)*(w*w+1-w*q/(1-q/2))
                    LC<-round(Mu+Lim*sqrt(Var*lambda/(2-lambda)),2)
                    
                    if(LC<=Min|LC<=Mu){
                              ARL0.ZS=1; ARL1.ZS=1;ARL0.SS=1; ARL1.SS=1 
                              salida<-c(ARL0.ZS,ARL1.ZS,ARL0.SS,ARL1.SS)
                              return(salida)}
                    
                    #####se construyen las matrices de probabilidades de transicion
                    m=max(1,ceiling((LC-Min)/amp))
                    mT=m*m #m estados por cada estadistico
                    
                    while (mT>m.max){
                              amp=amp+0.01       
                              m=ceiling((LC-Min)/amp)
                              mT=m*m        
                    }
                    
                    P=Q=matrix(0,ncol=mT,nrow=mT)
                    i.prom=floor(1+max((Mu-Min)/amp,0))
                    i.prom=(i.prom-1)*m+i.prom                    #Intervalo asociado a valor medio (estado bidimensional)
                    for(i.s in 1:m){                              # recorre todos los estados i para Ms
                              PMi.s=Min+(i.s-0.5)*amp                     # Punto medio del estado i.s
                              for(i.l in 1:m){                            # recorre todos los estados i para Ml
                                        PMi.l=Min+(i.l-0.5)*amp                   # Punto medio del estado i.l
                                        k=(i.l-1)*m+i.s                           ## transforma estados i.s,i.l a estado bidimensional k.
                                        for(h1 in 1:h){                           #Recorre todos el esapcio muestral, identificando el estado de llegada
                                                  Zs=lambda*Mst[h1]+(1-lambda)*PMi.s
                                                  Zl=lambda*Mlt[h1]+(1-lambda)*PMi.l
                                                  
                                                  if(max(Zs,Zl)-LC<(-5e-9)){#(-5e-9)
                                                            j.s=floor(1+max((Zs-Min)/amp,0))      #Obtiene el estado de destino j.s para cada Zts
                                                            j.l=floor(1+max((Zl-Min)/amp,0))      #Obtiene el estado de destino j.l para cada Ztl
                                                            l=(j.l-1)*m+j.s                       #transforma j.s, j.l a Estado bidimensional de llegada
                                                            if(l>mT){
                                                                      print(paste(LC,",",w,",",q,",",Lim,",",lambda))
                                                                      print(paste(Zs,",",Zl,",",Min,",",Mu,",",m,",",l,",",i.prom,",",nrow(P),ncol(P)))}
                                                            
                                                            P[k,l]=P[k,l]-Prob.IC[h1]             # Acumula las probabilidades en P[k,l]
                                                            Q[k,l]=Q[k,l]-Prob.OOC[h1]}           # Acumula las probabilidades en Q[k,l]
                                        }   #acumula en Negativo la probabilidades
                              } 
                    }
          }
          
          #Calculo del ARL - Cadena de Markov
          #if(opt==1){Id=diag(1,nrow=m)}else{Id=diag(1,nrow=mT)}
          if(opt==1){mT=m}else{mT=mT}
          
          for (i in 1:mT){
                    P[i,i]=1+P[i,i]
                    Q[i,i]=1+Q[i,i]    #Ahora P y Q contienes a Id-P e Id - Q respectivamente
          }
          if(image==T){
                    x11()
                    par(mfrow=c(1,2))
                    SparseM::image(as.matrix.csr(P))}
          
          if(rcond(P)<2.220446e-16){
                    ARL0.ZS=10000; ARL1.ZS=10000;ARL0.SS=10000; ARL1.SS=10000 
          }else{
                    P<-solve(P); Q=solve(Q)   # se reemplazan por las inversas, para evitar consumo de memoria
          
          if(image==T){SparseM::image(as.matrix.csr(P))}
          if(i.prom>mT){ARL0.ZS<-1; ARL1.ZS<-1;ARL0.SS<-1; ARL1.SS<-1;}
          else{
          ARL0.ZS=sum(P[i.prom,])
          ARL1.ZS=sum(Q[i.prom,])
          
          P.0=P[i.prom,]/ARL0.ZS
          
          ARL0.SS= sum(P.0%*%P)
          ARL1.SS= sum(P.0%*%Q)}} 
          salida<-c(ARL0.ZS,ARL1.ZS,ARL0.SS,ARL1.SS)
          return(salida)
}


#####Block 2 - Adaptive Univariate Gauge ##############################################################
#Functions for the evaluation of performance of proposed of based of gauge adaptive sample size CC    # 
#******************************************************************************************************

#ARL for adaptive Control Chart based on gauge
ARLwYsYl.TMV<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0){
  fact.n1<-factorial(n1)
  fact.n2<-factorial(n2)
  
  ARL=0
  LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;
  
  qd=conv(delta,r,q,mu,sig,dist,skew)[1:2]
  qs<-qd[1]; ql<-qd[2]
  
  p11=0; p12=0
  q11=0; q12=0
  
      for (s in 0:n1){
        fact.s=factorial(s)
        for (l in 0:(n1-s)){
          cte<-(fact.n1/(fact.s*factorial(l)*factorial(n1-s-l)))
          if (((w*s+l)<LW1)&((s+w*l)<LW1)){
            p11 = p11 + cte*(q/2)^(s+l)*((1-q)^(n1-s-l))
            q11 = q11 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))
            }
          else if(((w*s+l)<LC1)&((s+w*l)<LC1)){
            p12 = p12 + cte*(q/2)^(s+l)*((1-q)^(n1-s-l))
            q12 = q12 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))
            }
        }
      }
  
  p21=0; p22=0
  q21=0; q22=0
  
  for (s in 0:n2){
    fact.s=factorial(s)
      for (l in 0:(n2-s)){
        cte<-(fact.n2/(fact.s*factorial(l)*factorial(n2-s-l)))
        if (((w*s+l)<LW2)&((s+w*l)<LW2)){
          p21 = p21 + cte*(q/2)^(s+l)*((1-q)^(n2-s-l))
          q21 = q21 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))
          }
        else if(((w*s+l)<LC2)&((s+w*l)<LC2)){
        p22 = p22 + cte*(q/2)^(s+l)*((1-q)^(n2-s-l))
        q22 = q22 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))
        }
      }
    }

  p11=min(p11,0.999999999999999999999999999999);p12=min(p12,0.999999999999999999999999999999);
  p21=min(p21,0.999999999999999999999999999999);p22=min(p22,0.999999999999999999999999999999);
  
  q11=min(q11,0.999999999999999999999999999999);q12=min(q12,0.999999999999999999999999999999);
  q21=min(q21,0.999999999999999999999999999999);q22=min(q22,0.999999999999999999999999999999);
  
  p1=(1-p22)/(1-p22+p12)
  p2=p12/(1-p22+p12)
  
  #ARL[6]= (1-p22+p12)/((1-p11)*(1-p22)-p21*p12) ARL0 zero state
  ARL[1]= min(10000,(p1*(1-p22+p12)+p2*(1-p11+p21))/((1-p11)*(1-p22)-p21*p12)) #ARL0
  ARL[2]=  p1*n1+(1-p1)*n2    #ASS0
  ARL[3]= (p1*((1-p22)*n1+p12*n2)+p2*(p21*n1+(1-p11)*n2))/((1-p11)*(1-p22)-p21*p12) #Anos0
  ARL[4]= min(10000,(p1*(1-q22+q12)+p2*(1-q11+q21))/((1-q11)*(1-q22)-q21*q12)) #ARL1 
  ARL[5]= (p1*((1-q22)*n1+q12*n2)+p2*(q21*n1+(1-q11)*n2))/((1-q11)*(1-q22)-q21*q12) #ANos1
  
  names(ARL)<-c("ARL0 ", "E(n)0","ANOS0","ARL1","ANOS1" ) 
  return(round(ARL,3))
}

ARLwYsYl.DS<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0){
  
  #calcula las probabilidades en la colas, bajo y fuera de control
  qd=conv(delta,r,q,mu,sig,dist,skew)[1:2]
  qsd<-qd[1]; qld<-qd[2]
  
  
  #Probabilidad de no se?al en 2 es el complemento
  IC.PNS.n1<-0; OOC.PNS.n1<-0;
  IC.PS.n1<-0;  OOC.PS.n1<-0;
  IC.PS.n2<-0;  OOC.PS.n2<-0;
  n1.fact<-factorial(n1); n2.fact<-factorial(n2)
  
  for (s in 0:n1){
    s.fact<-factorial(s)
    for (l in 0:(n1-s)){
      wYSYL.n1<-w*s+l ; wYLYS.n1<-w*l+s;
      Cte1<-(n1.fact/(s.fact*factorial(l)*factorial(n1-s-l)))
      IC.Prob.n1 = Cte1*((q/2)^(s+l))*((1-q)^(n1-s-l))
      OOC.Prob.n1 = Cte1*(qsd^s)*(qld^l)*((1-qsd-qld)^(n1-s-l))
      
      if(wYSYL.n1<Lw&wYLYS.n1<Lw){IC.PNS.n1<-IC.PNS.n1+IC.Prob.n1;OOC.PNS.n1<-OOC.PNS.n1+OOC.Prob.n1}
      else if(wYSYL.n1>=Lc1|wYLYS.n1>=Lc1){IC.PS.n1<-IC.PS.n1+IC.Prob.n1;OOC.PS.n1<-OOC.PS.n1+OOC.Prob.n1}
      else{
        IC.Prob.n2<-0;OOC.Prob.n2<-0
        for (s2 in 0:n2){
          s2.fact<-factorial(s2)
          for (l2 in 0:(n2-s2)){
            wYSYL.n2<-w*s2+l2 ; wYLYS.n2<-w*l2+s2;
            Cte2<-(n2.fact/(s2.fact*factorial(l2)*factorial(n2-s2-l2)))
            if(wYSYL.n2>=Lc2-wYSYL.n1|wYLYS.n2>=Lc2-wYLYS.n1){
              IC.Prob.n2<-IC.Prob.n2+Cte2*((q/2)^(s2+l2))*((1-q)^(n2-s2-l2))
              OOC.Prob.n2<-OOC.Prob.n2 +Cte2*(qsd^s2)*(qld^l2)*((1-qsd-qld)^(n2-s2-l2))
            }
          }
        }
        IC.PS.n2 = IC.PS.n2 + IC.Prob.n2*IC.Prob.n1 
        OOC.PS.n2 = OOC.PS.n2 + OOC.Prob.n2*OOC.Prob.n1
      }
    }
  }
  
  ARL0<-1/(max(0.0001,IC.PS.n1+IC.PS.n2))
  ARL1<-1/(max(0.0001,OOC.PS.n1+OOC.PS.n2))
  ASS0<-n1+n2*(1-IC.PS.n1-IC.PNS.n1)
  ASS1<-n1+n2*(1-OOC.PS.n1-OOC.PNS.n1)
  ANOS0<-ARL0*ASS0
  ANOS1<-ARL1*ASS1 
  salida<-c(ARL0,ASS0,ANOS0,ARL1,ASS1,ANOS1)
  names(salida)<-c("ARL0","ASS0","ANOS0","ARL1","ASS1","ANOS1")
  return(round(salida,3))
}


#Profile ARl - Only Profile

CARLwYsYl.TMV<-function(n1,n2,Lc1,LW,Lc2,w,q,delta,r,curve="mean",mu=10,sig=1,dist="norm",skew=0){ 
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
    ARL[i]<-ARLwYsYl.TMV(n1,n2,Lc1,LW,Lc2,w,q,0,r[i],mu,sig,gamma,dist)[3]
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
    ARL[i]<-ARLwYsYl.TMV(n1,n2,Lc1,LW,Lc2,w,q,delta[i],1,mu,sig,gamma,dist)[3]
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
}


CARLwYsYl.DS<-function(n1,n2,Lc1,LW,Lc2,w,q,delta,r,curve="mean",mu=10,sig=1,dist="norm",skew=0){ 
  ARL=0
  
  if(curve=="dev"){
    for(i in 1:length(r)){
      ARL[i]<-ARLwYsYl.DS(n1,n2,Lc1,LW,Lc2,w,q,0,r[i],mu,sig,gamma,dist)[2]
    }
    return(cbind(r,ARL))
  }
  else if(curve=="mean"){
    for(i in 1:length(delta)){
      ARL[i]<-ARLwYsYl.DS(n1,n2,Lc1,LW,Lc2,w,q,delta[i],1,mu,sig,gamma,dist)[2]
    }
    return(cbind(delta,ARL))
  }
  
  else {stop("profile not enabled")}
}

#####Block 3 - Exactly optimizer for scheme based on gauge  ###########################################
#Optimizador de Graficos YT, YS-Yl, YS,YL, WYSYL y DSYSYL, con 1 y 3 puntos de cambio objetivos       #
#******************************************************************************************************

## Find q that meet ARL0 for(n,LC,w) configuration on the wYsYL Control Chart
RootARL.wYsYl<-function(n,Lc,w,ARLdes){
   Tol=0.001
   ninc<-0
   fact.n<-factorial(n)
   
   if (Lc>0.85*n){qi<-0.9}else if(Lc>0.5*n){qi<-0.7}else{qi<-(Lc/n)} ##ojo
   for (i in 1:2000){
      sumdf<-0;Prob<-0
      for(s in 0:n){
        fact.s<-factorial(s)
         for(l in 0:(n-s)){
            if((w*s+l<Lc)&(s+w*l<Lc)){
              cte<-(1/(fact.s*factorial(l)*factorial(n-s-l)))
               sumdf<-sumdf + cte*((qi/2)^(s+l-1))*((1-qi)^(n-s-l-1))*(s+l-n*qi)
               Prob = Prob + fact.n*cte*((qi/2)^s)*((qi/2)^l)*((1-qi)^(n-s-l))}
         }
      }
      fi<- (ARLdes-1)-ARLdes*Prob
      dfi<- (-ARLdes*fact.n/2)*sumdf
      q<- qi-fi/dfi
      if(q<0){q<-0.001; ninc<-ninc+1}else if(q>1){q<-0.99999; ninc<-ninc+1}
      if(abs(q-qi)<Tol|(ninc>3)){break};
      qi<-q      }
   return(q)                              }

## Find q that Meet ARL0, for(n,LCDif,LCSum) configuration on the joint Ys - YL; YS + YL Control Chart
RootARL.DSYsYl<-function(n,LcD,LcS,ARLdes){
   Tol=0.000001
   fact.n<-factorial(n)
   
   if (2*LcD>n){qi<-0.5}else{qi<-(LcD/n)}
   for (i in 1:2000){
      sumdf<-0;Prob<-0
      for(s in 0:n){
        fact.s<-factorial(s)
         for(l in 0:(n-s)){
            if((-s+l<LcD)&(s-l<LcD)&(s+l<LcS)){
               sumdf<-sumdf + (1/(fact.s*factorial(l)*factorial(n-s-l)))*((qi/2)^(s+l-1))*((1-qi)^(n-s-l-1))*(s+l-n*qi)
               Prob = Prob + (fact.n/(fact.s*factorial(l)*factorial(n-s-l)))*((qi/2)^s)*((qi/2)^l)*((1-qi)^(n-s-l))}
         }
      }
      fi<- (ARLdes-1)-ARLdes*Prob
      dfi<- (-ARLdes*factorial(n)/2)*sumdf
      q<- qi-fi/dfi
      if(q<0){q<-0.001}else if(q>1){q<-0.99999}
      if(abs(q-qi)<Tol)break;
      qi<-q          }
   return(q)}

# Set of values of WYSYL for fixed ( w,n).
Genera.LC<-function(w,n,decreasing = FALSE){
   Yl<-rep((0:n),n+1)
   Ys<-rep((0:n),each=n+1)
   Yt<-Yl+Ys
   Yl<-Yl[Yt<=n]
   Ys<-Ys[Yt<=n]
   LC<-unique(apply(data.frame(w*Ys+Yl,Ys+w*Yl),1,max))
   LC<-LC[LC>0]
   return(sort(LC,decreasing))
}

#Optimizer of wYsYl and specific schemes (esq="wYSYL","YSYL","YT","YS-YL") for a shift (delta,r) in the process

Opt.YsYl<-function(n,delta,r,ARLdes,esq,mu=10,sig=1,dist="norm",skew=0){ 
  t <- proc.time()
  par<-c(n,delta,r,ARLdes)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  print(par)
  print("_________________________________")
  ARL1.min=ARLdes; LC1=0;
  
  
  dw=0.1
  w.min=-2
  if (esq=="wYsYl"){# Optimice w, LC, q for scheme wYSYL
    w.p=c(1,0,-1,seq(w.min,-dw,dw),seq(dw,(1-dw),dw))
    for (i in 1:length(w.p)){
      L=Genera.LC(w.p[i],n)
      for(j in 1:length(L)) {
        q<-RootARL.wYsYl(n,L[j],w.p[i],ARLdes)
        if(q>0.998){break}
        else{
          qd<-conv(delta,r,q,mu,sig,dist,skew)[1:2]
          A0=ARLwYsYl(n,L[j],w.p[i],q/2,q/2);
          A1=ARLwYsYl(n,L[j],w.p[i],qd[1],qd[2])
          if((abs(A0-ARLdes)<=0.8)&(A1<(ARL1.min-0.05))){
            ARL1.min<-A1;LC1<-L[j];w<-w.p[i];qopt<-q}
        }
      }
      
    }
  }
  
  else if ((esq=="YT")|(esq=="YsYl")|(esq=="Ys-Yl")){ # Optimice LC y q for scheme YSYL (w=0) or YT (w=1) or YS-YL(W=-1)
    L<-1:n; w<-sum(c(-1,0,1)*(esq==c("Ys-Yl","YsYl","YT")))
    for (i in 1:length(L)){
      q<-RootARL.wYsYl(n,L[i],w,ARLdes)
      qd<-conv(delta,r,q,mu,sig,dist,skew)[1:2]
      A0=ARLwYsYl(n,L[i],w,q/2,q/2)
      A1=ARLwYsYl(n,L[i],w,qd[1],qd[2])
      if((abs(A0-ARLdes)<=0.8)&(A1<(ARL1.min-0.05))){
        ARL1.min<-A1;LC1<-L[i];qopt<-q }
    }
  }
  
  else{stop("esq must be: wYsYl, YsYl, Ys-Yl or YT")}
 
  sol<-c(n,round(LC1,2),round(w,3),round(qopt,5),round(ARL1.min,3))
  names(sol)<-c("n","Lc","w","q","ARL1")
  
  t <- proc.time()-t
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  print(paste("Esquema: ", esq))
  print(paste("ARL0=", round(ARLwYsYl(n,LC1,w,qopt/2,qopt/2),4)))
  print(paste("ARL1=", round(ARL1.min,4)))
  
  return(sol)
}


#Optimizer of esq="wYSYL", "DSYSYL" for three shifts (delta,r) (delta,0) (r,0) in the process

Opt.MDYsYl<-function(n,delta,r,ARLdes,esq,mu=10,sig=1,dist="norm",skew=0,pm=2,pmix=2,pd=1){ 
   t <- proc.time()
   par<-c(n,delta,r,ARLdes)
   names(par)<-c("n","delta","r","ARL.0 deseado")
   print(par)
   print("_________________________________")
   min=ARLdes*(pm+pmix+pd)
   LC1=0
   LC2=0

   dw=0.1
   w.min=-2
   if (esq=="wYsYl"){# Optimice w, LC, q for scheme wYSYL
     w.p=c(1,0,-1,seq(w.min,-dw,dw),seq(dw,(1-dw),dw))
     for (i in 1:length(w.p)){
       L=Genera.LC(w.p[i],n)
       for(j in 1:length(L)) {
         q<-RootARL.wYsYl(n,L[j],w.p[i],ARLdes)
         if(q>0.998){break}
         else{
           qd<-matrix(0,ncol=2,nrow=3)
           qd[1,]<-conv(delta,1,q,mu,sig,dist,skew)[1:2]
           qd[2,]<-conv(delta,r,q,mu,sig,dist,skew)[1:2]
           qd[3,]<-conv(0,r,q,mu,sig,dist,skew)[1:2]
           A0=ARLwYsYl(n,L[j],w.p[i],q/2,q/2);
           ARL1=c(ARLwYsYl(n,L[j],w.p[i],qd[1,1],qd[1,2]),ARLwYsYl(n,L[j],w.p[i],qd[2,1],qd[2,2]),ARLwYsYl(n,L[j],w.p[i],qd[3,1],qd[3,2]))
           A1= sum(c(pm,pmix,pd)*ARL1)
           if((abs(A0-ARLdes)<=0.8)&(A1<(min-0.05))){
             min<-A1;LC1<-L[j];w<-w.p[i];qopt<-q;ARL1.opt<-ARL1}
         }
       }
       
     }
     sol<-c(n,round(LC1,2),round(w,3),round(qopt,5))
     names(sol)<-c("n","Lc","w","q")
     t <- proc.time()-t
     print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
     print(paste("Esquema: ", esq))
     print(paste("FO=", round(min,4)))
     print(paste("ARL0=", round(ARLwYsYl(n,LC1,w,qopt/2,qopt/2),4)))
     print(paste("ARL(mu,sig)=", round(ARL1.opt[2],4)))
     print(paste("ARL(mu)=",round(ARL1.opt[1],4)))
     print(paste("ARL(sig)=", round(ARL1.opt[3],4)))
  }
 
 else if (esq=="DSYsYl"){
      L1<-1:n
      for (i in 1:length(L1)){
         for(j in 1:min(n+1,2*i)){
            q<-RootARL.DSYsYl(n,L1[i],j,ARLdes)
            qd<-matrix(0,ncol=2,nrow=3)
            qd[1,]<-conv(delta,1,q,mu,sig,dist,skew)[1:2]
            qd[2,]<-conv(delta,r,q,mu,sig,dist,skew)[1:2]
            qd[3,]<-conv(0,r,q,mu,sig,dist,skew)[1:2]
            A0=ARLDSYsYl(n,L1[i],j,q/2,q/2);
            ARL1=c(ARLDSYsYl(n,L1[i],j,qd[1,1],qd[1,2]),ARLDSYsYl(n,L1[i],j,qd[2,1],qd[2,2]),ARLDSYsYl(n,L1[i],j,qd[3,1],qd[3,2]))
            A1= sum(c(pm,pmix,pd)*ARL1)
            if((abs(A0-ARLdes)<=0.8)&(A1<(min-0.05))){
              min<-A1;LC1<-L1[i];LC2<-j;qopt<-q;ARL1.opt<-ARL1 
            }
         }
      }
      sol<-c(n,LC1,LC2,round(qopt,5))
      names(sol)<-c("n","LcD","LcS","q")
      t <- proc.time()-t
      print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
      print(paste("Esquema: ", esq))
      print(paste("FO=", round(min,4)))
      print(paste("ARL0=", round(ARLDSYsYl(n,LC1,LC2,qopt/2,qopt/2),4)))
      print(paste("ARL(mu,sig)=", round(ARL1.opt[2],4)))
      print(paste("ARL(mu)=",round(ARL1.opt[1],4)))
      print(paste("ARL(sig)=", round(ARL1.opt[3],4)))
   }
   else{stop("Esquema no Valido")}
   return(sol)
}

#####Block 4 - GA optimizer for scheme based on gauge  ################################################
#GA Optimizer for D2YSYL, WYSYL-TMV, WYSYL-DS                                                         #
#******************************************************************************************************

#### Decode, fitness and GA Functi?n for D2YSYL GA with one shift or three shift

decode.D2YsYl <- function(string,l,maxim) {# l: Longitudes de las 2 cadenas (LC,q), Min, Max: Minimos y maximos (LC,w,q)
   string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
   q <- q.min+(binary2decimal(string[(l[1]+1):(l[1]+l[2])])/maxim[2])*(q.max-q.min)
   Lc.min<-n*q/(1-q); Lc.max<-2*n*(1-q)/q
   Lc <- Lc.min+(binary2decimal(string[1:l[1]])/maxim[1])*(Lc.max-Lc.min)
   return(c(Lc,q))
}

fitn.D2YsYl<-function(string,l,maxim,mobj){
   sol<-decode.D2YsYl(string,l,maxim)
   Lc<-sol[1]; q<-sol[2]
   ARL0=ARLD2YsYl(n,Lc,q,q/2,q/2)
   
   if(mobj==F){
   qd = conv(delta,r,q,mu,sig,dist,skew)[1:2]
   FObj1 = ARLD2YsYl(n,Lc,q,qd[1],qd[2])
   if(ARL0==Inf|FObj1==Inf){fitn<-0}else{
   #fitn<-10000-w1*((ARL0<ARL0obj)*6*abs(ARL0-ARL0obj)+(ARL0>=ARL0obj)*abs(ARL0-ARL0obj))-w2*FObj1}
     fitn<-10000-w1*6*(ARL0<ARL0obj)*(ARL0obj-ARL0)-w2*(FObj1/ARL.ref-1)}
     }
   else{
   qdm<-conv(delta,1,q,mu,sig,dist,skew)[1:2]
   qdmix<-conv(delta,r,q,mu,sig,dist,skew)[1:2]
   qdd<-conv(0,r,q,mu,sig,dist,skew)[1:2]
   FObj1<-pm*ARLD2YsYl(n,Lc,q,qdm[1],qdm[2])+pmix*ARLD2YsYl(n,Lc,q,qdmix[1],qdmix[2])+pd*ARLD2YsYl(n,Lc,q,qdd[1],qdd[2])
   if(ARL0==Inf|FObj1==Inf){fitn<-0}else{
   fitn<-10000-w1*((ARL0<ARL0obj)*6*abs(ARL0-ARL0obj)+(ARL0>=ARL0obj)*abs(ARL0-ARL0obj))-w2*FObj1/(pm+pmix+pd)}
   }
   fitn    
}

GA.D2YsYl<-function(n,delta,r,ARL0obj,mu,sig,dist,skew,mobj=T){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  q.min<-0.0001;q.max<-0.65; prec.q<-0.0001
  prec.L<-0.01 
  l1<-length(decimal2binary(1/prec.L))#L
  l2<-9
  l<-c(l1,l2)
  maxim<-0
  for(i in 1:length(l)){
    maxim[i]<-binary2decimal(rep(1,l[i]))}
  
  GA.sol<-ga(type = "binary", nBits = l[1]+l[2], fitness = fitn.D2YsYl,l=l,maxim=maxim,mobj=mobj,
             popSize = 400, maxiter = 100, run = 50, pcrossover = 0.9, pmutation= 0.25) #.pendiente modificar tipo d cruce y seleccion
  sol<-decode.D2YsYl(GA.sol@solution[1,],l,maxim)
  Lc<-sol[1]; q<-sol[2]
  qdmix = conv(delta,r,q,mu,sig,dist,skew)[1:2]
  ARL<-c(ARLD2YsYl(n,Lc,q,q/2,q/2),ARLD2YsYl(n,Lc,q,qdmix[1],qdmix[2]))
  print(par)
  print("_________________________________")
  t <- proc.time()-t
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  print(paste("Esquema: ", "D2YsYl"))
  print(paste("ARL0=", round(ARL[1],4)))
  print(paste("ARL(mu=",mu+delta*r,"sig=",sig*r,")=", round(ARL[2],4)))
  if(mobj==T){
    qdm<-conv(delta,1,q,mu,sig,dist,skew)[1:2]
    qdd<-conv(0,r,q,mu,sig,dist,skew)[1:2]
    ARL<-c(ARL,ARLD2YsYl(n,Lc,q,qdm[1],qdm[2]),ARLD2YsYl(n,Lc,q,qdd[1],qdd[2]))
    print(paste("ARL(mu=",mu+delta*r,"sig=1)=", round(ARL[3],4)))
    print(paste("ARL(mu=0, sig=",sig*r,")=", round(ARL[4],4)))
    sol<-c(n,sol,ARL)
    names(sol)<-c("n","LC","q","ARL0","ARL1.mix","ARL1.mu","ARL1.dev")
  }
  else{sol<-c(n,sol,ARL)
  names(sol)<-c("n","LC","q","ARL0","ARL1")}
  
  return(sol)
}

#### Decode, fitness and GA Functi?n for wYSYL TMV with one shift

decode.wYsYlTMV <- function(string,n2.max,maxim,l,F.Lc1,F.w,F.n1) {#l:longitud de las cadenas,
   #string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
   
   q <- q.min + (q.max-q.min)*(binary2decimal(string[1:l[1]]))/maxim[1]   #(max(q.min,q.max-binary2decimal(string[1:l[1]])*prec.q)
   if(is.na(F.w)==F){w=F.w}else{
   w <- w.min + (w.max-w.min)*(binary2decimal(string[(l[1]+1):sum(l[1:2])]))/maxim[2]  #max(w.min,w.max-(binary2decimal(string[(l[1]+1):sum(l[1:2])]))*prec.w)
   }
   
   if(is.na(F.n1)==F){n1=F.n1}else{
             n1 <- 1+ floor((binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])])/maxim[3])*(n-2))
             #n1 <- max(ceiling((n-1)*binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])])/maxim[3]),1)# max(1, n-1-binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])]))  #esto le da mas probabilidad a los numeros peque?os
   }
   #n2 <- n + max(ceiling((n2.max-n)*binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])/maxim[4]),1)#max(n+1,n2.max-binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])) #esto le da mas probabilidad a los numeros peque?os
   n2 <- (n+1)+ floor((binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])/maxim[4])*(n.max-n-1))
   Lw <- 0.05 + (1-0.05)*(binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])]))/maxim[5]    #max(LC.min, 1-binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])])*prec.LC)
   Lc.min<- Lw +0.05
   if(is.na(F.Lc1)==F){Lc1=F.Lc1}else{
   Lc1<-Lc.min + (1.05-Lc.min)*binary2decimal(string[(sum(l[1:5])+1):sum(l[1:6])])/maxim[6]        #max(Lc.min, 1.05 - binary2decimal(string[(sum(l[1:5])+1):(sum(l[1:5])+long)])*prec.LC)
   }
   Lc2<-Lc.min + (1-Lc.min)*binary2decimal(string[(sum(l[1:6])+1):sum(l)])/maxim[7]      #min(1, Lc.min+binary2decimal(string[(sum(l[1:6])+1):(sum(l))])*prec.LC)
   return(c(q,w,n1,n2,Lw,Lc1,Lc2))
}

fitn.wYsYlTMV<-function(string,n2.max,maxim,l,F.Lc1,F.w,F.n1,fit){
   sol<-decode.wYsYlTMV(string,n2.max,maxim,l,F.Lc1,F.w,F.n1)
   q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
   ARL=ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
   if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
   else{
   if(fit==1){ fitn<-10000-w1*((ARL[1]<ARL0obj)*(ARL0obj-ARL[1]))-w3*(ARL[2]>n)*(ARL[2]-n)-w2*(ARL[4]/ARL.ref-1)}
     else{fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref-1)}
    }
   fitn
}

GA.wYsYlTMV<-function(n,delta,r,ARL0obj,n2.max=2*n,mu,sig,dist,skew,maxiter,F.Lc1=NA,F.w=NA,F.n1=NA,fit=2){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  prec.w<-0.05
  prec.q<-0.0001
  prec.L<-0.01 ###Fundamental para la precision del algoritmo.....(0.05, mejor opcion)
  
  if(delta==0){l1<-9}else{l1<-10}
  
  #l1<-length(decimal2binary((q.max-q.min)/prec.q))#q
  
  if(is.na(F.w)==F){l2<-0}else{
  l2<-length(decimal2binary((w.max-w.min)/prec.w +1))#w
    #l2<-5
  }
  if(is.na(F.n1)==F){l3<-0}else{
            l3<-length(decimal2binary(n-2))#n1
            #l3<-length(decimal2binary(n-1+1))#n1
  }
  
  l4<-length(decimal2binary(n2.max-n-1))#n2
  
  l5<-length(decimal2binary(1.05/prec.L))#Lw
  if(is.na(F.Lc1)==F){l6<-0}else{
  l6<-length(decimal2binary(1/prec.L))#LC1
  }
  l7<-length(decimal2binary(1/prec.L))#LC2
  l<-c(l1,l2,l3,l4,l5,l6,l7)
  
  maxim<-0
  for(i in 1:length(l)){
  maxim[i]<-binary2decimal(rep(1,l[i]))}#calcula los maximos que se obtienen con cada longitud de gen
  
  
  GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlTMV,l=l,maxim=maxim,F.Lc1=F.Lc1,F.w=F.w,n.max,F.n1,fit=fit,
             popSize = 500, maxiter = maxiter, run = 100, pcrossover = 0.95, pmutation= 0.25) #.pendiente modificar tipo d cruce y seleccion
  
  sol<-decode.wYsYlTMV(GA.sol@solution[1,],n.max,maxim,l,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
  sol.ARL=ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,sol.ARL),4)
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q","ARL0 ", "E(n)0","ANOS0","ARL1","ANOS1" )
  return(salida)
}

#### Decode, fitness and GA Functi?n for wYSYL DS with one shift

decode.wYsYlDS <- function(string,n.max,maxim,l,F.Lc1,F.w,F.n1) {#l:longitud de las cadenas,
  string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
  
  q <- q.min + (q.max-q.min)*binary2decimal(string[1:l[1]])/maxim[1]   #(max(q.min,q.max-binary2decimal(string[1:l[1]])*prec.q)
  if(is.na(F.w)==F){w=F.w}else{
    w <- w.min + (w.max-w.min)*binary2decimal(string[(l[1]+1):sum(l[1:2])])/maxim[2]  #max(w.min,w.max-(binary2decimal(string[(l[1]+1):sum(l[1:2])]))*prec.w)
  }
  
  if(is.na(F.n1)==F){n1=F.n1}else{
   n1<-1+ floor((binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])])/maxim[3])*(n-2))
  #n1 <- max(ceiling((n-1)*binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])])/maxim[3]),1)     # max(1, n-1-binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])]))  #esto le da mas probabilidad a los numeros peque?os
  }
  n2<-(n-n1)+ floor((binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])])/maxim[3])*(n.max-n))
  #n2 <- (n-n1)+max(ceiling((n.max-n)*binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])/maxim[4]),1)#max(n+1,n2.max-binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])) #esto le da mas probabilidad a los numeros peque?os
  
  Lw <- 0.01+(n1-0.01)*binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])])/maxim[5]    #max(LC.min, 1-binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])])*prec.LC)
  Lc.min<- Lw +0.05
  if(is.na(F.Lc1)==F){Lc1=n1+0.5}else{
    Lc1<-Lc.min + (n1+0.5-Lc.min)*binary2decimal(string[(sum(l[1:5])+1):sum(l[1:6])])/maxim[6]        #max(Lc.min, 1.05 - binary2decimal(string[(sum(l[1:5])+1):(sum(l[1:5])+long)])*prec.LC)
  }
  Lc2<-Lc.min +(n1+n2-Lc.min)*binary2decimal(string[(sum(l[1:6])+1):sum(l)])/maxim[7]      #min(1, Lc.min+binary2decimal(string[(sum(l[1:6])+1):(sum(l))])*prec.LC)
  return(c(q,w,n1,n2,Lw,Lc1,Lc2))
}

fitn.wYsYlDS<-function(string,n.max,maxim,l,F.Lc1,F.w,F.n1,fit){
  sol<-decode.wYsYlDS(string,n.max,maxim,l,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
  ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  else{
    if(fit==1){fitn<-10000-w1*max(ARL0obj-ARL[1],0)-w3*max(ARL[2]-n,0)-w2*(ARL[4]/ARL.ref-1)}
    else{fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref-1)}
  }
  fitn
  
}

GA.wYsYlDS<-function(n,delta,r,ARL0obj,n.max=2*n,mu,sig,dist,skew,maxiter,F.Lc1=NA,F.w=NA,F.n1=NA,fit=2){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  
  w.min<-(-2);w.max<-1.1; prec.w<-0.05
  q.min<-0.001;q.max<-0.65; prec.q<-0.005
  prec.L<-0.1 
  
  if(delta==0){l1<-9}else{l1<-10}
  
  #l1<-length(decimal2binary((q.max-q.min)/prec.q))#q

  if(is.na(F.w)==F){l2<-0}else{
    l2<-length(decimal2binary((w.max-w.min)/prec.w+1))#w
    #l2<-5
  }
  if(is.na(F.n1)==F){l3<-0}else{
    l3<-length(decimal2binary(n-1+1))#n1
  }
  
  l4<-length(decimal2binary(n.max-n+1))#n2
  
  #l5<-length(decimal2binary(1/prec.L))#Lc1
  l5<-length(decimal2binary(min(10*(n-1)+1,101)))
  if(is.na(F.Lc1)==F){l6<-0}else{
    #l6<-length(decimal2binary(1/prec.L))#Lw
    l6<-length(decimal2binary(min(10*(n-1)+1,101)))#Lw
  }
  #l7<-length(decimal2binary(1/prec.L))#LC2
  l7<-length(decimal2binary(min(100+1,1/prec.L+1)))#LC2
  l<-c(l1,l2,l3,l4,l5,l6,l7)
  
  maxim<-0
  for(i in 1:length(l)){
    maxim[i]<-binary2decimal(rep(1,l[i]))}#calcula los maximos que se obtienen con cada longitud de gen
  
  
  GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlDS,l=l,maxim=maxim,F.Lc1=F.Lc1,F.w=F.w,n.max,F.n1=F.n1,fit=fit,
             popSize = 500, maxiter = maxiter, run = 100, pcrossover =0.95 , pmutation= 0.25) #.pendiente modificar tipo d cruce y seleccion
  
  sol<-decode.wYsYlDS(GA.sol@solution[1,],n.max,maxim,l,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
  sol.ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,sol.ARL),5)
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q","ARL0","E(n)0","ANOS0","ARL1","E(n)1","ANOS1")
  return(salida)
}

#### Decode, fitness and GA Functión for wYSYL EWMA with one shift

decode.wYsYlEWMA <- function(string,maxim,l,F.w){#l:longitud de las cadenas,
          string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
          
          q <- round(q.min + (q.max-q.min)*binary2decimal(string[1:l[1]])/maxim[1],5)   #(max(q.min,q.max-binary2decimal(string[1:l[1]])*prec.q)
          
          if(is.na(F.w)==F){w=F.w}else{w <- round(w.min + (w.max-w.min)*binary2decimal(string[(l[1]+1):sum(l[1:2])])/maxim[2],2)}
          Lam<-0.02+floor(binary2decimal((string[(sum(l[1:2])+1):sum(l[1:3])]))*49/maxim[3])*0.02
          #Lam<-round(0.05+(0.95)*binary2decimal((string[(sum(l[1:2])+1):sum(l[1:3])]))/maxim[3],2)
          Lim<-round(0.05+(3.45)*binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])/maxim[4],4)
          
          return(c(q,w,Lam,Lim))
}

fitn.wYsYlEWMA<-function(string,maxim,l,F.w,fit,ARL.ref){
          sol<-decode.wYsYlEWMA(string,maxim,l,F.w)
          q<-sol[1];w<-sol[2];Lam<-sol[3];Lim<-sol[4]
          ARL=ARLwYsYl.EWMA(n,Lim,w,Lam,q,delta,r,opt,mu,sig,dist,skew)
          if(ARL[3]==Inf|ARL[4]==Inf|ARL[3]<1|ARL[4]<1){fitn<-0}
          else{
                    if(fit==1){fitn<-10000-w1*max(ARL0obj-ARL[3],0)-w2*(ARL[4]/ARL.ref-1)}
                    else{fitn<-10000-w1*((ARL[3]<ARL0obj)*2*(ARL0obj-ARL[3])+(ARL[3]>=ARL0obj)*(ARL[3]-ARL0obj))-w2*(ARL[4]/ARL.ref-1)}
          }
          fitn
}


GA.wYsYlEWMA<-function(n,delta,r,ARL0obj,mu,sig,dist,skew,maxiter,F.w=NA,fit=2){
          t <- proc.time() 
          param<-c(n,delta,r,ARL0obj)
          names(param)<-c("n","delta","r","ARL.0 deseado")
          
          prec.w<-0.05 ; prec.q<-0.0001; prec.Lim<-0.01; prec.Lam<-0.05 
          if(delta==0){l1<-10}else{l1<-11}
          #l1<-length(decimal2binary((q.max-q.min)/prec.q))#q
          if(is.na(F.w)==F){l2<-0}else{l2<-length(decimal2binary((w.max-w.min)/prec.w))}#w
          l3<-length(decimal2binary(49))
          l4<-length(decimal2binary(3.5/prec.Lim))
          
          l<-c(l1,l2,l3,l4)
          
          maxim<-0
          for(i in 1:length(l)){maxim[i]<-binary2decimal(rep(1,l[i]))}#calcula los maximos que se obtienen con cada longitud de gen
          
          ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
          
          GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlEWMA,l=l,maxim=maxim,F.w=F.w,fit=fit,ARL.ref=ARL.ref,
                     popSize = 250, maxiter = maxiter, run = 20, pcrossover =0.95 , pmutation= 0.25,parallel = F) #.pendiente modificar tipo d cruce y seleccion
          
          #sol<-decode.wYsYlEWMA(GA.sol@solution[1,],maxim,l,F.w)
          #q<-sol[1];w<-sol[2];Lam<-sol[3];Lim<-sol[4]
          #ARL.ref=ARLwYsYl.EWMA(n,Lim,w,Lam,q,delta,r,opt,mu,sig,dist,skew)[4]
          #GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlEWMA,l=l,maxim=maxim,F.w=F.w,fit=fit,ARL.ref=ARL.ref,
                    # popSize = 300, maxiter = maxiter,suggestions=GA.sol@solution[1,], run = 20, pcrossover =0.95 , pmutation= 0.25,parallel = F) #.pendiente modificar tipo d cruce y seleccion
          
          sol<-decode.wYsYlEWMA(GA.sol@solution[1,],maxim,l,F.w)
          q<-sol[1];w<-sol[2];Lam<-sol[3];Lim<-sol[4]
          sol.ARL=ARLwYsYl.EWMA(n,Lim,w,Lam,q,delta,r,opt,mu,sig,dist,skew)
          salida<-c(n,Lim,w,Lam,q,sol.ARL[3:4])
          t <- proc.time()-t
          print(param)
          print("_________________________________")
          print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
          names(salida)<-c("n","Lim","w","Lambda","q","ARL0","ARL1")
          return(salida)
}


#####Block 5 - Simulator for scheme based on gauge     ################################################
#Simulator for wYSYL - WYSYL-TMV and wYSYL-DS                                                         #
#******************************************************************************************************

#wYsYl Simulator 
SARLwYsYl<-function(n,Lc,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,it=10000){
  
 param<-conv(delta,r,q,mu,sig,dist,skew)

  Xi= param[3];Xs= param[4]

  if(dist=="norm"){mu1=mu+delta*sig; sig1=r*sig}
  else if(dist=="snorm"){xi1=param[8];omega1=param[9];lambda1=param[10]}
  else{gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
  LR=0
  
  for (j in 1:it){
    EC="True"; rl=0;  
    
    while (EC=="True"){
      rl=rl+1;  contyl=0; contys=0;
      if(dist=="norm"){Dat= rnorm(n,mu1, sig1)}
      else if(dist=="snorm"){Dat= rsn(n,xi1,omega1,lambda1)}
      else{Dat= rlnorm3(n,gamma1,mulog1,siglog1)}
      
      contys= sum(Dat>Xs); contyl= sum(Dat<Xi);
      
      if((contyl+w*contys < Lc) & (w*contyl+ contys< Lc)){
        (EC="True")}else (EC="False")
    }
    LR[j]=rl;
  }
  arl=sum(LR)/length(LR)
  sdrl=sd(LR)
  sol<-c(round(arl,3),round(sdrl,3))
  if(delta==0&r==1){et<-0}else{et<-1}
  names(sol)<-paste(c("ARL","SDRL"),et)
  return(sol)
}

#wYSYL.TMV Simulator 

SARLwYsYl.TMV<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,it=10000){#Simulador wYSYL.TMV
  
  LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;

  param<-conv(delta,r,q,mu,sig,dist,skew)
  Xi= param[3];Xs= param[4]
  
  if(dist=="norm"){mu1=mu+delta*sig; sig1=r*sig}
  else if(dist=="snorm"){xi1=param[8];omega1=param[9];lambda1=param[10]}
  else{gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
   LR=0;Np=0;
   
   p11= Fmultw(n1,LW1,w,q/2, q/2)
   p12= Fmultw(n1,LC1,w,q/2, q/2)-p11
   p21= Fmultw(n2,LW2,w,q/2, q/2)
   p22= Fmultw(n2,LC2,w,q/2, q/2)-p21
   
   p1=(1-p22)/(1-p22+p12)
   p2=p12/(1-p22+p12)
   
   for (j in 1:it){
      EC="True"; rl=0;  ncum=0; f2=0;
      if (runif(1)<= p1) (f2=0) else (f2=1);

      while (EC=="True"){
         n =(1-f2)*n1+f2*n2;
         rl=rl+1; ncum=ncum+n; contyl=0; contys=0;
         LC=(1-f2)*LC1 + f2*LC2; LA=(1-f2)*LW1 + f2*LW2;
         if(dist=="norm"){Dat= rnorm(n,mu1, sig1)}
         else if(dist=="snorm"){Dat= rsn(n,xi1,omega1,lambda1)}
         else{Dat= rlnorm3(n,gamma1,mulog1,siglog1)}
         contys= sum(Dat>Xs); contyl= sum(Dat<Xi);
         
         if((contyl+w*contys < LC) & (w*contyl+ contys< LC)){
            (EC="True")
            if((contyl+w*contys < LA) & (w*contyl+ contys< LA)) (f2=0) else (f2=1)
         }
         else (EC="False")
      }
      LR[j]=rl; Np[j]=ncum; #
   }
   arl=sum(LR)/length(LR)
   Nprom= sum(Np/LR)/length(Np)
   sdarl=sd(LR)
   ATI= sum(Np)/length(Np)
   sol<-c(round(arl,3),round(Nprom,3),round(sdarl,3),round(ATI,3))
   if(delta==0&r==1){et<-0}else{et<-1}
   names(sol)<-paste(c("ARL", "E(n)","SDRL","ATI"),et)
   return(sol)
}

#wYSYL-DS Simulator

SARLwYsYl.DS<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,it=10000){#Simulador wYsYl
  
  param<-conv(delta,r,q,mu,sig,dist,skew)
  
  S= param[3];L= param[4]
  
  if(dist=="norm"){mu1=mu+delta*sig; sig1=r*sig}
  else if(dist=="snorm"){xi1=param[8];omega1=param[9];lambda1=param[10]}
  else{gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
  LR=0;Np=0;
  
  
  for (j in 1:it){
    EC="True"; rl=0; ncum=0
    while (EC=="True"){
      rl=rl+1;
      if(dist=="norm"){Dat= rnorm(n1+n2,mu1, sig1)}
      else if(dist=="snorm"){Dat= rsn(n1+n2,xi1,omega1,lambda1)}
      else{Dat= rlnorm3(n1+n2,gamma1,mulog1,siglog1)}
      
      ncum=ncum+n1
      YL1= sum(Dat[1:n1]>L); YS1= sum(Dat[1:n1]<S);
      
      if((YL1+w*YS1 < Lw) & (w*YL1+ YS1< Lw)){EC="True"}
      else if((YL1+w*YS1 < Lc1) & (w*YL1+ YS1< Lc1)){
        ncum = ncum+n2
        YL2= sum(Dat>L); YS2= sum(Dat<S);
        if((YL2+w*YS2 < Lc2) & (w*YL2+ YS2< Lc2)){EC="True"}else{EC="False"}
        }
      else {EC="False"}
    }
    LR[j]=rl;Np[j]=ncum;
  }
  arl=sum(LR)/length(LR)
  Nprom= sum(Np/LR)/length(Np)
  sdarl=sd(LR)
  ANOS= sum(Np)/length(Np)
  sol<-c(round(arl,3),round(Nprom,3),round(sdarl,3),round(ANOS,3))
  if(delta==0&r==1){et<-0}else{et<-1}
  names(sol)<-paste(c("ARL", "E(n)","SDRL","ANOS"),et)
  return(sol)
}

#wYsYl- EWMA Simulator
SARLwYsYl.EWMA<-function(n,Lim,w,lambda,q,delta,r,opt=1,mu=10,sig=1,dist="norm",skew=0,it=30000){#Simulador wYsYl
  
  #Obtiene los valores de S y L
  param<-conv(delta,r,q,mu,sig,dist,skew)
  S= param[3];L= param[4]
  
  #Encuentra los parametros de la distribución fuera de control del proceso
  if(dist=="norm"){
    mu1=mu+delta*sig; sig1=r*sig}else if(dist=="snorm"){
      xi1=param[8];omega1=param[9];lambda1=param[10]}else{
        gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
  ##ojoooooo, no esta incluid la distribución weibull
  
  # Obtiene el valor esperado para Z0
  fact.n=factorial(n)
  
  if(opt==1){
    Mu=0; E2=0   #Primer y segundo momento no central
    for (s in 0:n){
      fact.s<-factorial(s)
      for (l in 0:(n-s)){
        p<-fact.n/(fact.s*factorial(l)*factorial(n-s-l))*(q/2)^(s+l)*((1-q)^(n-s-l))
        Mt=max(w*s+l,s+w*l)
        Mu=Mu+Mt*p        #Valor esperado
        E2=E2+(Mt^2)*p    #Esperanza del cuadrado
      }}
    Var=E2-(Mu^2)          #varianza
  }else{Mu=n*(q/2)*(w+1);Var=n*q/2*(1-q/2)*(w*w+1-w*q/(1-q/2))}    # Media y valor esperado para opción 2.
  
  LC<-round(Mu+Lim*sqrt(Var*lambda/(2-lambda)),2)     #Factor L se convierte en Limite de Control LC
  
  ##Simula el RL
  LR=0                                   #incia vector para almacenar las it longitudes de racha
  for (j in 1:it){                        # it iteraciones
    EC="True"; rl=0;                      #proceso inicia bajo control, con longitud de racha 0
    Zt=Mu; Zlt=Mu; Zst=Mu                 # el EWMA empieza en su valor central
    while (EC=="True"){
      rl=rl+1
      if(dist=="norm"){Dat= rnorm(n,mu1,sig1)}else if(dist=="snorm"){
        Dat= rsn(n,xi1,omega1,lambda1)}else{
          Dat= rlnorm3(n,gamma1,mulog1,siglog1)}     ##ojoooooo, no esta incluid la distribución weibull
      
      Yl= sum(Dat>L); Ys= sum(Dat<S);           #Obtiene los conteos
      Mst=w*Ys+Yl; Mlt=Ys+w*Yl; Mt=max(Mst,Mlt) # Calcula los estadisticos
      
      if(opt==1){
        Zt=lambda*Mt+(1-lambda)*Zt
      }else{
        Zlt=lambda*Mlt+(1-lambda)*Zlt
        Zst=lambda*Mst+(1-lambda)*Zst
        Zt=max(Zlt,Zst)}
      if(Zt-LC<(-5e-9)){EC="True"}else{EC="False"}
    }
    LR[j]=rl                                    # Almacena valor de la longitud de racha e inicia nueva iteración
  }
  
  arl=sum(LR)/length(LR)          
  sdrl=sd(LR)
  sol<-c(round(arl,3),round(sdrl,3))            # Calcula El ARL y SDRL
  if(delta==0&r==1){et<-0}else{et<-1}
  names(sol)<-paste(c("ARL", "SDRL"),et)
  return(sol)
}


####################################### End #############################################################