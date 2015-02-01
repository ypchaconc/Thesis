
# Graphics and Chains

# Hyperameters

pdf(file="/home/mirt/Documentos/Thesis/R/Chains/Hyperparameters.pdf",width=12,height=8.5)
par(mfrow=c(2,3))

plot(as.ts(salida$alpha.theta.smpl),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$alpha.theta.smpl),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$alpha.theta.smpl),cex.main=0.9, col="red", lwd=2)


plot(as.ts(salida$beta.theta.smpl),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$beta.theta.smpl),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$beta.theta.smpl),cex.main=0.9, col="red", lwd=2)

dev.off()


# Latent Trait

x11()
par(mfrow=c(1,3))

plot(as.ts(salida$theta.smpl[,100]),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$theta.smpl[,1]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$theta.smpl[,1]),cex.main=0.9, col="red", lwd=2)



postscript(file="/home/mirt/Documentos/Thesis/R/Chains/Triats.eps", 
           paper = "special", family = "ComputerModern", encoding = "TeXtext.enc", 
           onefile=FALSE, horizontal=FALSE,width=8,height=12)

# Items parameters

nitems<-dim(salida$alpha.smpl)[2]
nind<-dim(salida$theta)[2]
pdf(file="/home/mirt/Documentos/Thesis/R/Chains/Items.pdf",width=12,height=8.5)
par(mfrow=c(2,3))
for(i in 1:nitems){
  plot(as.ts(salida$alpha.smpl[,i]),cex.main=0.9,ylab="",xlab="iterations")
  autocorr.plot(mcmc(salida$alpha.smpl[,i]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
  plot(density(salida$alpha.smpl[,i]),cex.main=0.9, col="red", lwd=2)

  plot(as.ts(salida$beta.smpl[,i]),cex.main=0.9,ylab="",xlab="iterations")
  autocorr.plot(mcmc(salida$beta.smpl[,i]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
  plot(density(salida$beta.smpl[,i]),cex.main=0.9, col="red", lwd=2)  
}
dev.off()

