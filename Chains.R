
# Graphics and Chains

# Hyperameters

x11()
par(mfrow=c(2,3))

plot(as.ts(salida$alpha.theta.smpl),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$alpha.theta.smpl),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$alpha.theta.smpl),cex.main=0.9, col="red", lwd=2)


plot(as.ts(salida$beta.theta.smpl),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$beta.theta.smpl),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$beta.theta.smpl),cex.main=0.9, col="red", lwd=2)

# Latent Trait

x11()
par(mfrow=c(1,3))

plot(as.ts(salida$theta.smpl[,1]),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$theta.smpl[,1]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$theta.smpl[,1]),cex.main=0.9, col="red", lwd=2)

# Items parameters

x11()
par(mfrow=c(2,3))

plot(as.ts(salida$alpha.smpl[,1]),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$alpha.smpl[,1]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$alpha.smpl[,1]),cex.main=0.9, col="red", lwd=2)


plot(as.ts(salida$beta.smpl[,1]),cex.main=0.9,ylab="",xlab="iterations")
autocorr.plot(mcmc(salida$beta.smpl[,1]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=1000)
plot(density(salida$beta.smpl[,1]),cex.main=0.9, col="red", lwd=2)
