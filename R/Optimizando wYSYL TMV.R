
######                  WYsYl.TMV                      ######     
######        Test  - Junio 2018                       ######

source('~/Doctorado-EyO/Tesis Doctorado/Script R/GaugeSPCFunctions.R', encoding = 'UTF-8')


### Optimización esquemas wYsYl (FSS-VSS-DS-EWMA)

# Parametros del proceso bajo control
mu<-10; sig<-1; skew<-0; dist<-"norm"

#Parametros diseño Grafico de Control
n<-7; delta<-0.6; r=1; ARL0obj=370

#Parametro Adicionales para VSS - DS -EWMA.
n.max<-2*n; opt=1; m.max=1000
ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
#Parametors para algoritmo genetico.
w1<-1 ; w2 <-200; w3<-100 # pesos de las 3 componentes de la funcion fitness 
pm<-2;pmix<-3;pd<-1

q.min<-0.01
if(delta==0){q.max<-0.30}else{q.max=0.65}
w.min<-(-2);w.max<-1

###Optimización
Opt.YsYl(n,delta,r,ARL0obj,"YT",mu,sig,dist,skew)
Opt.YsYl(n,delta,r,ARL0obj,"Ys-Yl",mu,sig,dist,skew)
Opt.YsYl(n,delta,r,ARL0obj,"YsYl",mu,sig,dist,skew)
Opt.YsYl(n,delta,r,ARL0obj,"wYsYl",mu,sig,dist,skew)
Opt.MDYsYl(n,delta,r,ARL0obj,"wYsYl",mu,sig,dist,skew,pm,pmix,pd)
Opt.MDYsYl(n,delta,r,ARL0obj,"DSYsYl",mu,sig,dist,skew,pm,pmix,pd)
GA.D2YsYl(n,delta,r,ARL0obj,mu,sig,dist,skew,mobj=F)
#Opt.wYSYL.TMV(n,delta,r,ARL0obj,n.max,mu,sig,dist,skew)
#Opt.wYSYL.DS(n,delta,r,ARL0obj,F.w=NA,n.max,mu,sig,dist,skew)
GA.wYsYlTMV(n,delta,r,ARL0obj,mu,sig,dist,skew,maxiter=800,n2.max=n.max,F.Lc1=1.05,F.w=-1,F.n1=NA,fit=2)
GA.wYsYlDS(n,delta,r,ARL0obj,mu,sig,dist,skew,maxiter=800,n.max=n.max,F.Lc1=T,F.w=NA,F.n1=NA,fit=2)
GA.wYsYlEWMA(n,delta,r,ARL0obj,mu,sig,dist,skew,100,F.w=-1,fit=2)


############Verificación de un plan , simulación vs calculo###########

it<-50000

####Comparacion wYSYL
n<-10 ; LC<-4 ; w=-1; q=0.5710;
c(ARLwYsYl(n,LC,w,q/2,q/2),ARLwYsYl(n,LC,w,qd[1],qd[2]))
c(SARLwYsYl(n,LC,w,q,0,1,mu,sig,dist,skew,it)[1],SARLwYsYl(n,LC,w,q,delta,r,mu,sig,dist,skew,it)[1])

####Comparacion wYSYL.TMV- stady state

n1=5; n2=15; Lc1=1.09; Lw=0.3; Lc2=0.5; w=-1; q=0.5710
ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
c(SARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,0,1,mu,sig,dist,skew,it)[1:2],SARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,it)[c(1,4)])

####Comparacion wYSYL.DS -steady state
n0=5;
n1=3; n2=10; Lw=2; Lc1=2; Lc2=7.41;w=-0.35; 
q=Root.wYSYL.DS(w, n1, n2, n0, Lw, Lc1, tol=0.001,maxiter=100)[1]
ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
c(SARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,0,1,mu,sig,dist,skew,it)[1:2],SARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,it))

####Comparacion wYSYL.EWMA - Zero state
Lim=2.878;w=-1;lambda=0.17;q=0.36302

ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,mu,sig,dist,skew,image=T)
c(SARLwYsYl.EWMA(n,Lim,w,lambda,q,0,1,opt,mu,sig,dist,skew,it)[1],SARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,mu,sig,dist,skew,it)[1])



####Comparacion wYSYL.EWMA - Zero state
n=3; Lim=0.312 ;w=0.71;lambda=0.05;q=0.00309;delta=0.2;r=1;mu=10;sig=1;dist="norm";skew=0;opt=1;
it=10
ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,mu,sig,dist,skew,image=F)
c(SARLwYsYl.EWMA(n,Lim,w,lambda,q,0,1,opt,mu,sig,dist,skew,it)[1],SARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,mu,sig,dist,skew,it)[1])











##############################################################################
# Prueba Algortimo Búsqueda Total wYsYl.TMV - Febrero 2017
##############################################################################
h<-function(ARL0,w,n0, n1, n2, LW1, LW2, LC1, LC2,q){
  fact.n1<-factorial(n1)
  fact.n2<-factorial(n2)          
  qs<-q/2; ql<-q/2
  p11=0; p12=0
  
  for (s in 0:n1){
    fact.s<-factorial(s)
    for (l in 0:(n1-s)){
      cte<-(fact.n1/(fact.s*factorial(l)*factorial(n1-s-l)))
      if (((w*s+l)<LW1)&((s+w*l)<LW1)){p11 = p11 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))}
      else if(((w*s+l)<LC1)&((s+w*l)<LC1)){p12 = p12 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))}
    }
  }
  
  p21=0; p22=0
  for (s in 0:n2){
    fact.s<-factorial(s)
    for (l in 0:(n2-s)){
      cte<-(fact.n2/(fact.s*factorial(l)*factorial(n2-s-l)))
      if (((w*s+l)<LW2)&((s+w*l)<LW2)){p21 = p21 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))}
      else if(((w*s+l)<LC2)&((s+w*l)<LC2)){p22 = p22 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))}
    }
  }
  
  
  p1= (1-p22)/(1-p22+p12)
  n=p1*n1+(1-p1)*n2
  ARL=(1-p11+p21+p1*(p11+p12-p21-p22))/((1-p11)*(1-p22)-p21*p12)
  f.p1= ((n/n0)-1)^2
  f.ARL<-((ARL/ARL0)-1)^2
  f.glob<-f.p1+f.ARL
  return(c(f.glob,f.p1,f.ARL))
}

Root.wYSYL.TMV<-function(ARL0, w, n0, n1, n2, LW1, LW2, LC1, LC2,tol=0.001,maxiter=100){
  
  R=(sqrt(5)-1)/2; qmin=0.01;qmax=0.7  #Golde search, busca el minimo de la funcion de perdida (restricciones)
  k=0
  
  
  repeat{
    k=k+1
    hmax=h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,qmax)[1]
    if(hmax<0.7|k==7){qmax=qmax+0.1;break}else{qmax=qmax-0.1}
    
  }
  
  
  q1=qmax-R*(qmax-qmin);q2=qmin +R*(qmax-qmin)
  h1=h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,q1)[1]
  h2=h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,q2)[1]
  hmin=h(ARL0, w,n0, n1, n2, LW1, LW2, LC1, LC2,qmin)[1]
  
  if((hmin<h1)&(hmax<h2)){qopt=NA;hopt=NA}
  else{
    repeat{
      k=k+1
      
      if(h1<h2){
        hopt=h1
        qopt=q1
        qmax=q2
        q2=q1;h2=h1
        q1=qmax-R*(qmax-qmin)
        h1=h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,q1)[1]
        
      }
      else{
        hopt=h2
        qopt=q2
        qmin=q1 
        q1=q2;h1=h2
        q2=qmin +R*(qmax-qmin)
        h2=h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,q2)[1]
      }
      
      if ((R*(qmax-qmin)/qopt)<=tol|k>=maxiter){break}
    }} 
  
  return(c(qopt,hopt,k,(R*(qmax-qmin)/qopt)))
}


Opt.wYSYL.TMV<-function(n0,delta,r,ARLdes,n.max=4*n0,mu=10,sig=1,dist="norm",skew=0){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  
prec.w<-0.1; w.min<-(-2.0)
Omega.w<-c(1,-1,0,seq(1-prec.w,prec.w,-prec.w),seq(-prec.w,-1+prec.w,-prec.w),seq(-1-prec.w,w.min,-prec.w))            #k
ARL1.min=ARLdes
delta.ARL1<-2 ; delta.n<-0.1
cr<-0  #cuenta el número de raices halladas, soluciones factibles.

for(w in Omega.w){
for (n1 in (1:(n0-1))){
  Omega.Lw1<-(Genera.LC(w,n1)-0.001)/n1
  for(n2 in ((n0+1):n.max)){
    Omega.Lw2<-(Genera.LC(w,n2)-0.001)/n2
     L<-sort(unique(c(Omega.Lw1,Omega.Lw2)))
     Omega.Lw<-L[L<1]
     for(Lw in Omega.Lw){#no puede lw=1,1.05
       Omega.Lc1<- c(Omega.Lw[Omega.Lw>Lw],1,1.05)
        for(Lc1 in Omega.Lc1){##desde el proximo valor a lw hasta 1.05
          Omega.Lc2<-Omega.Lc1[-length(Omega.Lc1)]
          for(Lc2 in Omega.Lc2){##desde el proximo valor a lw hasta 1. 1.05 solo para n1
            Root<-Root.wYSYL.TMV(ARLdes,w,n0, n1, n2, Lw*n1, Lw*n2, Lc1*n1, Lc2*n2)
            if(is.na(Root[1])==F&&Root[2]<=(0.05^2)){
              sol<-ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,Root[1],delta,r) 
              if((sol[1]> (ARLdes-delta.ARL1))&(sol[2]<=(n0+delta.n))){
              cr<-cr+1 #Solucion factible
              # print(c(cr,Root[1],w,n1,n2,Lc1,Lw,Lc2,sol[1:3])) Imprime la solucion factible 
              if(sol[4]<ARL1.min){
                ARL1.min<-sol[4]
                Opt<-c(n1,n2,Lc1,Lw,Lc2,w,Root[1],sol)
                }}
              }
            }
          }
        }
      }
    }
  }

t <- proc.time()-t
print(par)
print("_________________________________")
print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))

Opt[8:12]<-ARLwYsYl.TMV(Opt[1],Opt[2],Opt[3],Opt[4],Opt[5],Opt[6],Opt[7],delta,r,mu,sig,dist,skew)

salida<-c(Opt,cr)
names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q","ARL0","E(n)0","ANOS0","ARL1","ANOS1","#sol")
return(salida)

}


#####Grafico prueba para verificar que la raiz es hallada
q<-seq(0.01,6,0.005)

prueb=0
for(i in 1:length(q)){
  prueb[i]<- h(ARL0,w, n0, n1, n2, LW1, LW2, LC1, LC2,q[i])[1]
}

x11()
plot(q,prueb,type="l",xlim=c(0.01,1),ylim=c(0,2))
abline(h=0)
abline(v=2.997874e-01)
q[prueb==min(prueb)]


###################################################################
#####grafico de la de distribución de wYSYL 
###################################################################


w=0.33
n=30
q=0.3
delta=0.5; r=1
qs=ql=q/2

x11()
#○par(mfrow=c(2,1))

Dist<-pwYSYL(w,n,q/2,q/2)
plot(Dist[,1],Dist[,2],type="h", xlab=expression(w*Y[S]*Y[L]),ylab="Probabilidad")
text(0.8*n, 0.8*max(Dist[,2]),paste("(n = ",n, "; w = ",w,"q = ",qs+ql,")"),cex=1.3)

qd=convnorm(delta,r,q)
Dist<-pwYSYL(w,n,qd[1],qd[2])
points(Dist[,1],Dist[,2],type="h", col="red")


##########################Validación óptimos con genético########################
##se alimenta con tabla de parametros xls


###verificando arl0 y nprom, se toman 12 columnas

param<-read.table("clipboard",header=F,sep="\t",dec=",")

res<-matrix(0,ncol=6,nrow=nrow(param))

for (i in 1:nrow(param)){
  n1<-param[i,11]; n2<-param[i,12]
  Lc1<-param[i,7];Lc2<-param[i,8];LW<-param[i,9];
  w<-param[i,10];q<-param[i,6]
  delta<-param[i,5];r<-param[i,4];
  res[i,]<-ARLwYsYl.DS(n1,n2,Lc1,LW,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm")
}

res
write.table(res,"res.txt")

#######################probando óptimo entre filas 

res<-matrix(0,ncol=nrow(param),nrow=nrow(param))

for (i in 1:nrow(param)){
  n1<-param[i,11]; n2<-param[i,12]
  Lc1<-param[i,7];Lc2<-param[i,8];LW<-param[i,9];
  w<-param[i,10];q<-param[i,6]
  for(j in 1:nrow(param)){
    delta<-param[j,5];r<-param[j,4];
    res[j,i]<-ARLwYsYl.DS(n1,n2,Lc1,LW,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm")[4]
  }
}

ind<-1:ncol(res)
res<-cbind(res,rep(0,nrow(res)))
for(i in 1:nrow(res)){
  res[i,ncol(res)]<-min(ind[min(res[i,-ncol(res)])==res[i,-ncol(res)]])
}


res


#########################Calculando los ARL1 optimos para wYSYL, para contraste

param<-read.table("clipboard",header=F,sep="\t",dec=",")

ARL1<-0
for(i in 1:nrow(param)){
n<-param[i,3];ARLdes<-param[i,2]
delta<-param[i,5];r<-param[i,4];
ARL1[i] <-Opt.YSYL(n,delta,r,ARLdes,"wYSYL",mu=10,sig=1,gam=0,dist="norm")[5]
}

param<-cbind(param,ARL1)

write.table(param,"ARL.txt")


####dibujando ARL1 para los esquema óptimos 

n<-10; ARLdes<-370
param<-read.table("clipboard",header=F,sep="\t",dec=",")

delta<-seq(0,2.5,0.05)
r<-seq(1,2.5,0.02)
curv="mean"

x11()
par(mfrow=c(1,2))

tit<-paste("n1","n2","Lc1","LW","Lc2","q","w",collapse=",",sep=" ; ")

leg<-"a"

parwYsYl<-Opt.YsYl(n,delta[2],r[2],ARLdes,"wYsYl",mu=10,sig=1,gam=0,dist="norm")[2:4]

for (i in 1:nrow(param)){#Escala Original ARL
  q<-param[i,1]
  Lc1<-param[i,2];Lc2<-param[i,3];LW<-param[i,4];
  w<-param[i,5];
  n1<-param[i,6]; n2<-param[i,7]
  CARL<-CARLwYsYl.TMV(n1,n2,Lc1,LW,Lc2,w,q,delta,r,curve=curv,mu=10,sig=1,gamma=0,dist="norm")
  
  leg[i]<-paste(n1,n2,Lc1,LW,Lc2,q,w,collapse=",",sep=" ; ")
  
  
  
  if(i==1){   
    plot(CARL[,1],CARL[,2],type="l",axes=F,xlab="",ylab="",ylim=c(0,max(CARL[,2])),lwd=2)  
    axis(1,at=seq(min(CARL[,1]),max(CARL[,1]),0.2),labels=seq(min(CARL[,1]),max(CARL[,1]),0.2),pos=0,cex.axis=1.5)
    axis(2,c((0:5)*100),c((0:5)*100), pos=min(CARL[,1]),cex.axis=1.5)
    mtext("ARL",2,line=1.5,cex=1.5)
    if (curv=="mean"){mtext(expression(delta),1,line=1,cex=1.5)}else{mtext(expression(italic("r")),1,line=1.5,cex=1.5)}
     }
  else{
    lines(CARL[,1],CARL[,2],lty=i,lwd=2) 
  }
}
 legend(min(CARL[,1])+0.4,200,leg,bty="n",lty=1:nrow(param),lwd=2,title=tit,title.adj=0.25,cex=1)
 CARL<-CARLwYsYl(n,parwYsYl[1],parwYsYl[2],parwYsYl[3],delta,r,curve=curv,mu=10,sig=1,gamma=0,dist="norm")
 lines(CARL[,1],CARL[,2],lty=1,lwd=2, col="red") 

for (i in 1:nrow(param)){   #Escala Logaritmica
q<-param[i,1]
Lc1<-param[i,2];Lc2<-param[i,3];LW<-param[i,4];
w<-param[i,5];
n1<-param[i,6]; n2<-param[i,7]
CARL<-CARLwYsYl.TMV(n1,n2,Lc1,LW,Lc2,w,q,delta,r,curve=curv,mu=10,sig=1,gamma=0,dist="norm")

if(i==1){
  plot(CARL[,1],log(CARL[,2]),type="l",axes=F,xlab="",ylab="",ylim=c(0,log(max(CARL[,2]))),lwd=2)  
  axis(1,at=seq(min(CARL[,1]),max(CARL[,1]),0.2),labels=seq(min(CARL[,1]),max(CARL[,1]),0.2),pos=0,cex.axis=1.5)
  axis(2,-1:10,-1:10, pos=min(CARL[,1]),cex.axis=1.5)
  mtext("Log(ARL)",2,line=1.5,cex=1.5)
  if (curv=="mean"){mtext(expression(delta),1,line=1,cex=1.5)}else{mtext(expression(italic("r")),1,line=1.5,cex=1.5)}
}
else{
lines(CARL[,1],log(CARL[,2]),lty=i,lwd=2) 
}

}
 
legend(min(CARL[,1])+0.4,200*log(max(CARL[,2]))/max(CARL[,2]),leg,bty="n",lty=1:nrow(param),lwd=2,title=tit,title.adj=0.25,cex=1)
CARL<-CARLwYsYl(n,parwYsYl[1],parwYsYl[2],parwYsYl[3],delta,r,curve=curv,mu=10,sig=1,gamma=0,dist="norm")
lines(CARL[,1],log(CARL[,2]),lty=1,lwd=2, col="red") 




###########grafico de control wYSYL - DS  --- esquema general

CL1=5.5; CL2=7.5 ; WL=3
X1=c(2, 3, 7, 5)
X2=c(NA,6,NA,9)

#jpeg("Figure1.jpeg")
x11()
plot(-1,-1,ylim=c(0,15),xlim=c(0.5,4),yaxt="n",xaxt="n",ylab="", xlab=expression("Time"),cex.axis=1.2,cex.lab=1.4)
points(1:length(X1),X1,pch=c(20,1,15,1),cex=1.7)
lines(1:length(X1),X1,lw=2)
points(c(2,4),X2[c(2,4)],pch=c(20,15),cex=1.7)
abline(h=c(WL,CL1,CL2),lty=c(2,1,1),lw=3)
lines(c(2,2),c(X1[2],X2[2]),lty=2,lwd=2)
lines(c(4,4),c(X1[4],X2[4]),lty=2,lwd=2)
text(0.65,CL1+0.5,expression(italic("UCL")^"I"),cex=1.3)
text(0.65,CL2+0.5,expression(italic("UCL")^"II"),cex=1.3)
text(0.65,WL+0.5,expression(italic("WL")),cex=1.3)
mtext(expression(italic("Max")(italic("wY")[S]+italic("Y")[L],italic("Y")[S]+italic("wY")[L])),2,line=1.5,cex=1.4,las=0)
axis(1,at=1:4,labels=c("A","B","C","D"),cex=1.4)

legend(2.5,14,c("In control"," 2nd Stage Required","Out of control"),title=expression(bold("Process Diagnosis")),pch=c(20,1,15),bty="n")

###########grafico de control wYSYL - DS  --- ejemplo


n1=5; n2= 7 ; w=0; q= 0.1630; CL1=5.5; CL2=5 ; WL=2
mu<-4
sigma<-0.3; delta=1.0 ; r=1.4

k = -qnorm(q/2)
S=mu-k*sigma; L=mu+k*sigma
q1= c(pnorm((-k-delta)/r),1-pnorm((-k+delta)/r))

#X0<-matrix(rnorm(10*(n1+n2),mu,sigma),nrow=n1+n2,ncol=10)
#X1<-matrix(rnorm(1*(n1+n2),mu+delta*sigma,sigma*r),nrow=n1+n2,ncol=1)
#X<-cbind(X0,X1)
#write.table(X,"Xsample.txt") guardado en carpeta del articulo

setwd("~/Doctorado-EyO/Tesis Doctorado/2. Grafico wYSYL-TMV/Articulo -versiones")
X<-read.table("Xsample.txt")
X[,4]<-X[,1]
X[,1]<-X[,3]
X[,7]<-X[,2]
YS1<-0  ;YL1<-0
YST<-0  ;YLT<-0
WYSYL1<-0;  WYSYLT<-0


pc=c(4,13,19)
pc=c(20,1,15)

#jpeg("Figure3.jpeg")
x11()
plot(-1,-1,ylim=c(0,n1+n2),xlim=c(1,ncol(X)),yaxt="n",ylab="", xlab=expression("Time"),cex.axis=1.2,cex.lab=1.4)
abline(h=c(WL,CL2),lty=c(2,1),lw=3)

for(j in 1:ncol(X)){
  
  YS1[j]<-sum(X[1:n1,j]<S)   ;YL1[j]<-sum(X[1:n1,j]>L)
  WYSYL1[j]<-max(w*YS1[j]+YL1[j],YS1[j]+w*YL1[j])
  pc1=pc[1]
  if(WYSYL1[j]<WL){WYSYLT[j]<-NA}
  else{
    pc1=pc[2]
    pc2=pc[1]
    YST[j]<-sum(X[,j]<S)   ;YLT[j]<-sum(X[,j]>L)
    WYSYLT[j]<-max(w*YST[j]+YLT[j],YST[j]+w*YLT[j])
    if(WYSYLT[j]>=CL2){pc2=pc[3]}
    points(j,WYSYLT[j],pch=pc2,cex=1.7)
    lines(c(j,j),c(WYSYL1[j],WYSYLT[j]),lty=2,lwd=2)
   }
    points(j,WYSYL1[j],pch=pc1,cex=1.7)
 }
  lines(1:ncol(X),WYSYL1,lw=2)


  text(2,CL2+0.5,expression(italic("UCL")^"II" == 5),cex=1.3)
  text(2,WL+0.5,expression(italic("WL") == 2),cex=1.3)
  mtext(expression(italic("Max")(italic("Y")[S],italic("Y")[L])),2,line=1.5,cex=1.7,las=0)

legend(7.5,11,c("In control"," 2nd Stage Required","Out of control"),title=expression(bold("Process Diagnosis")),pch=pc,bty="n")

dev.off()


#########################EWMA WYSYL ############################################






ARLs<-c(SARLwYsYl.EWMA(n,LC,w,lambda,q,0,1,opt,mu,sig,dist,skew,it=30000)[1],SARLwYsYl.EWMA(n,LC,w,lambda,q,delta,r,opt,mu,sig,dist,skew,it=30000)[1])

mp=c(10,20,30,50,75,100,150,200,250,300,400,450,500,550,600,800,1000,2000,4000,6000,8000)
ARLc<-matrix(0,ncol=2,nrow=length(mp))

i=0
for(m in mp){
  i=i+1
  if(opt==2){madj=ceiling(sqrt(m))}else{madj=m}
  ARLc[i,]<-ARLwYsYl.EWMA(n,LC,w,lambda,q,delta,r,madj,opt,mu,sig,dist,skew)[1:2]
}

x11()
plot(mp,ARLc[,1],type="l",lty=1,ylim=c(0,400))
abline(h=ARLs,col=c("Black","Gray"),lty=2)
lines(mp,ARLc[,2],col="Gray",lty=2)











#####Steiner comparacion
n=1
mu=0; sig=1; J=0.5
w=-1; q=pnorm(-1)*2; lambda=0.2045; Lim=2.78
dist="norm"

delta=J/sqrt(n); r=1;opt=1

(ARL<-ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt))


#####################

# print everything after loop is finished
for (i in 0:101) {
          print(i)
          Sys.sleep(0.01)
}

# simplist way to print within loop
for (i in 0:101) {
          print(i)
          Sys.sleep(0.01)
          flush.console()
}

# fancier text progress
# install.packages("svMisc")
require(svMisc)
for (i in 0:101) {
          progress(i)
          Sys.sleep(0.01)
          if (i == 101) cat("Done!\n")
}

# fancier text progress with bar
# install.packages("svMisc")
require(svMisc)
for (i in 0:101) {
          progress(i, progress.bar = TRUE)
          Sys.sleep(0.01)
          if (i == 101) cat("Done!\n")
}



#####Prueba EWMA sobre 50 casos aleatorios.
n.casos=50
m.n<-sample(c(5,10,15,20),n.casos,replace=T)
m.Lim<-round(runif(n.casos,0.05,3.2),3)
m.w<-round(runif(n.casos,-2,1),2)
m.lambda<-round(runif(n.casos,0.05,1),2)
m.q<-round(runif(n.casos,0.05,0.50),5)
m.delta<-round(runif(n.casos,0.5, 1.5),2)
m.r<-round(runif(n.casos,1.2, 2),2)

res<-matrix(0, nrow=n.casos,ncol=17)
colnames(res)<-c("n","Lim","w","Lambda","q","delta","r","Time.O1.C","A0.O1C","A1.O1C","A0.O1S","A1.O1S","Time.O2.C","A0.O2C","A1.O2C","A0.O2S","A1.O2S")

for(i in 1:n.casos){
          
          n=m.n[i]; Lim=m.Lim[i];w=m.w[i];lambda=m.lambda[i];q=m.q[i]; delta=m.delta[i]; r=m.r[i];
          
          
          res[i,1:7] <-c(n,Lim,w,lambda,q,delta,r)
          opt=1
          res[i,8] <-system.time(ARL<-ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt))[3]
          res[i,9:10] <-ARL[1:2] 
          res[i,11:12] <-c(SARLwYsYl.EWMA(n,Lim,w,lambda,q,0,1,opt,it=50000)[1],SARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,it=10000)[1]) 
          opt=2
          res[i,13] <-system.time(ARL<-ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt))[3]
          res[i,14:15] <-ARL[1:2] 
          res[i,16:17] <-c(SARLwYsYl.EWMA(n,Lim,w,lambda,q,0,1,opt,it=50000)[1],SARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,opt,it=10000)[1]) 
}


###Comparacion casos opcin 2 con m=4000

casos<- read.table("clipboard", header=T,dec=",")
opt=2
res<-matrix(0, nrow=nrow(casos),ncol=3)
colnames(res)<-c("Time.O2.C","A0.O2S","A1.O2S")



for(i in 1:nrow(casos)){
          
          n=casos$n[i]; Lim=casos$L[i];w=casos$w[i];lambda=casos$l[i];q=casos$q[i]; delta=casos$delta[i]; r=casos$r[i];
          res[i,1] <-system.time(ARL<-ARLwYsYl.EWMA(n,Lim,w,lambda,q,delta,r,m,opt))[3]
          res[i,2:3] <-ARL[1:2] 
}




# Parametros del proceso bajo control
mu<-10; sig<-1; skew<-0; dist<-"norm"

#Parametors para algoritmo genetico.
w1<-1 
w.min<-(-2);w.max<-1

casos<-read.table("clipboard", header=T,dec=",")
opt=1;w2=200;
res<-matrix(0,nrow=nrow(casos),ncol=10)
colnames(res)<-c("n","Lim","w","Lambda","q","ARLo.opt","ARL1.opt","ARL1.alt","tiempo.opt","tiempo.alt")

for(i in 1:nrow(casos)){
          n=casos$n0[i]; delta=casos$d[i]; r=casos$r[i]; ARL0obj=casos$ARL0[i];cambio=casos$Tcambio[i]
          #ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
          ARL.ref=casos$ARL.WYSYL.DS[i]
          
          if(n<=3){q.min<-0.002}else{q.min<-0.001}
          if(delta==0){q.max<-0.30}else{q.max=0.65}
          time<-system.time(sol1<-GA.wYsYlEWMA(n,delta,r,ARL0obj,mu,sig,dist,skew,50,F.w=NA,fit=2))[3]
          
          if(cambio=="Media"){Fw=-1}else if(cambio=="Varianza"){Fw=1}else{Fw=0}
          timef<-system.time(sol2<-GA.wYsYlEWMA(n,delta,r,ARL0obj,mu,sig,dist,skew,50,F.w=Fw,fit=2))[3]
          
          if(sol1[7]<sol2[7]){sol<-sol1}else{sol<-sol2}
          res[i,1:7] <-sol
          res[i,8] <-max(sol1[7],sol2[7]) 
          res[i,9:10]<-c(time,timef)
          write.table(res,"res.txt")
}
library(doParallel)
library(parallel)
library(foreach)
library(iterators)
