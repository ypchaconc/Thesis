rm(list = ls())
setwd("/home/mirt/Documentos/Thesis/R/Simulaciones")
library("VGAM")
x<-seq(0,1,length.out=100)###to plot
n<-15 ##numero ind
nitem<-8 ###numero items
ntest<-10 ####numero de tests
wt<-runif(1,.2,.8)
at<-runif(1,0,5)
bt<-log(.5)/log(1-wt^(at))
plot(x,dkumar(x,at,bt),type="l")
abline(v=qkumar(.15,at,bt))
abline(v=qkumar(.85,at,bt))

theta<- rkumar(n, shape1= at,shape2= bt)#simular trazo
summary(theta)
#################
###items
#################
(w<-runif(nitem,qkumar(.15,at,bt),qkumar(.85,at,bt)))
alpha<-runif(nitem,0,10)
beta<-log(.5)/log(1-w^(alpha))
(pitem<-cbind(alpha,beta))
##(pitem<-matrix(rexp(nitem*2,1),ncol=2))##simular parametros de item
################
###CCI
################
plot(x,pkumar(x,alpha[1],beta[1]),type="l",col=rainbow(length(beta))[1])
for(i in 2:length(beta)){
  points(x,pkumar(x,alpha[i],beta[i]),type="l",col=rainbow(length(beta))[i])
}

####matrices para calcular P(Yij=1|theta y pitem)
tm<-theta%*%t(rep(1,nitem))
pima<-rep(1,n)%*%t(pitem[,1])
pimb<-rep(1,n)%*%t(pitem[,2])
####P(Yij=1|theta y pitem)
ps<-apply(matrix(1:(n*nitem),ncol=nitem),1:2,function(i)pkumar(tm[i],pima[i],pimb[i]))
rm(list = c("tm","pima","pimb"))####removiendo las matrices temporales
###para generar los test y tenerlos disponibles en la lista rtest
rtest<-lapply(1:ntest,function(j)ifelse(matrix(runif(n*nitem),ncol=nitem)<ps,1,0))
####quitando individuos con todo 0 o 1 se recomienda no usar
for(i in 1:length(rtest)){
  cnt<-apply(rtest[[i]],1,function(x)sum(x)/dim(rtest[[i]])[2])
  names(cnt)<-1:dim(rtest[[i]])[1]
  rtest[[i]]<-rtest[[i]][-c(as.numeric(names(cnt[cnt==0 | cnt==1]))),]
}
###para generar los test y exportarlos a archivo de texto
a<-lapply(1:ntest,function(j)write.table(ifelse(matrix(runif(n*nitem),ncol=nitem)<ps,1,0), 
                                         file = paste(j,"_kstest.txt",sep=""), sep = "\t"))
rm("a")###removiendo a
###para comprobar, si se compilo rtest

mean(sapply(1:ntest,function(i)rtest[[i]][8,4]))
ps[8,4]

