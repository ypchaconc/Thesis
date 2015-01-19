

# Generador datasets

rm(list = ls())
setwd("/home/mirt/Documentos/Thesis/R/Simulaciones")

library("VGAM")


# numero ind
n<-1000 
# numero items
nitem<-50
# numero de tests
ntest<-2

set.seed(200)
wt<-runif(1, .2, .8)

at<-runif(1, 1, 10)
# bt<-runif(1, 1, 10)
bt<-log(.5)/log(1-wt^(at))

# to plot
x<-seq(0,1,length.out=100)
plot(x,dkumar(x,at,bt),type="l")
abline(v=qkumar(.15,at,bt))
abline(v=qkumar(.85,at,bt))

# simular trazo
theta<- rkumar(n, shape1= at,shape2= bt)
print("theta")
print(theta)
summary(theta)

plot(density(theta))
(meanKumar<- (bt *(gamma(1 + 1/at))*gamma(bt))/(gamma(1 + bt + 1/at) ) )

#################
# items
#################

(w<-runif(nitem,qkumar(.15,at,bt),qkumar(.85,at,bt)))
alpha<-runif(nitem, 1, 10)
# beta<- runif(nitem, 1, 10)
beta<-log(.5)/log(1-w^(alpha))
(pitem<-cbind(alpha,beta))
##(pitem<-matrix(rexp(nitem*2,1),ncol=2))##simular parametros de item


################
###CCI
################
# Poner labels a estas curvas para analizar comportamiento

plot(x,pkumar(x,alpha[1],beta[1]),type="l",col=rainbow(length(beta))[1])
for(i in 2:length(beta)){
  points(x,pkumar(x,alpha[i],beta[i]),type="l",col=rainbow(length(beta))[i])
}

x11()
plot(x,dkumar(x,alpha[1],beta[1]),type="l",col=rainbow(length(beta))[1], ylim = c(0, 15), )
for(i in 2:length(beta)){
  points(x,dkumar(x,alpha[i],beta[i]),type="l",col=rainbow(length(beta))[i])
 }


####matrices para calcular P(Yij=1|theta y pitem)
tm<-theta%*%t(rep(1,nitem))
pima<-rep(1,n)%*%t(pitem[,1])
pimb<-rep(1,n)%*%t(pitem[,2])

####P(Yij=1|theta y pitem)
ps<-apply(matrix(1:(n*nitem),ncol=nitem),1:2,function(i)pkumar(tm[i],pima[i],pimb[i]))

####removiendo las matrices temporales
rm(list = c("tm","pima","pimb"))

###para generar los test y tenerlos disponibles en la lista rtest
rtest<-lapply(1:ntest,function(j)ifelse(matrix(runif(n*nitem),ncol=nitem)<ps,1,0))

####quitando individuos con todo 0 o 1 se recomienda no usar
for(i in 1:length(rtest)){
  cnt<-apply(rtest[[i]],1,function(x)sum(x)/dim(rtest[[i]])[2])
  names(cnt)<-1:dim(rtest[[i]])[1]
  rtest[[i]]<-rtest[[i]][-c(as.numeric(names(cnt[cnt==0 | cnt==1]))),]
}

###para generar los test y exportarlos a archivo de texto
for(i in 1:ntest){
  
  #genera Tabla
  tabla = ifelse(matrix(runif(n*nitem),ncol=nitem)<ps,1,0)
  
  #scores nulos y perfectos
  inds = which(rowSums(tabla) == 0)
  if(!length(inds) == 0){
    minAlpha = min(alpha)
    posMinAlpha = which(alpha == minAlpha)
    tabla[inds,posMinAlpha] = 1
  }
  
  inds = which(rowSums(tabla) == ncol(tabla))
  if(!length(inds) == 0){
    maxAlpha = max(alpha)
    posMaxAlpha = which(alpha == maxAlpha)
    tabla[inds,posMaxAlpha] = 0
  }
  
  #exporta tabla
  write.table(tabla, 
              file = paste(i,"_kstest.csv",sep=""), sep = " ",row.names=F)
}



###para comprobar, si se compilo rtest

mean(sapply(1:ntest,function(i)rtest[[i]][8,4]))
ps[8,4]


sink("/home/mirt/Documentos/Thesis/R/param.txt", append = F)
at
bt
pitem
theta
sink()

