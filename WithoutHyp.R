###################################
# Universidad Nacional de Colombia
# Yuli Chacon
# Director A. Montenegro
# Algoritmo MH
# Modelo KK
###################################

# utilities libraries


library("VGAM")
library("truncnorm")
library(car)
library(coda)


rm(list = ls())

x<-read.table(file="/home/mirt/Documentos/Thesis/R/Simulaciones/1_kstest.txt" ,header=T,sep="")
x<-as.matrix(x)

# sink("/home/mirt/Documentos/Thesis/R/Simulaciones/Salida.txt")

KKModel = function (x, mcmc_size=10, show.iteration=TRUE,
                   
                    tao.theta = 1, 
                    tao.alpha = 0.25, 
                    tao.beta=0.25, 
                    
                    var.alpha = 100, 
                    var.beta = 100){
  
  # x: test data
  # mcmc_size: number of samples
  # file.name: file name to save the results, recomended)
  # show.iteration=TRUE: show each iteration
  # lambda.alpha.theta: prior distribution parameter hyperparemeter alpha.theta
  # lambda.beta.theta: prior distribution parameter hyperparemeter beta.theta 
  # alpha.theta: prior distribution parameter theta
  # beta.theta: prior distribution parameter theta
  # tao.theta: d  tunning parameter for latent traits (proposal)
  # tao.alpha: tunning parameter for the alpha parameter
  # tao.beta: tunning parameter for the beta parameter
  # lambda.alpha: parameter prior distribution
  # lambda.beta: parameter prior distribution
  
  
  #############################################################
  #               control section
  #############################################################
  if (missing(x))
  {stop(" array x containig the sample must be specified...")}
  
  ##################################################################
  #              global constants
  ##################################################################
  # N: sample size
  # I: test size
  N <- dim(x)[1]
  I <- dim(x)[2]
  
  set.seed(100)
  ##################################################################
  #           global variables
  ##################################################################
  
  # arrays to save the samples with inital values
  alpha.smpl <- matrix(1,nrow = mcmc_size, ncol = I)
  beta.smpl <- matrix(0,nrow = mcmc_size, ncol = I)
  theta.smpl <- matrix(0,nrow = mcmc_size, ncol = N)
  
  # Other global values
  
  # aceptation rate for thetas
  accept.rate.t <- rep(0, N)
  # aceptation rate for items
  accept.rate.i <- 0
 
  ##################################################################
  #           initial values
  ##################################################################
  # First values in the chains
  
  # alpha.theta.smpl[1] <-  rgamma(1, shape = 1, rate = lambda.alpha.theta)+1
  alpha.theta.smpl <-  4.938939
  
  # beta.theta.smpl[1]  <-  rgamma(1, shape = 1, rate = lambda.beta.theta)+1
  beta.theta.smpl <- 6.696828
  
  # theta.smpl[1,] <-  rkumar(N, 2, 2.5)
  theta.smpl[1,] <- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/theta.csv")))
  
  
  # alpha.smpl[1,] <-   rgamma(I, shape =1 , rate = lambda.alpha)
  alpha.smpl[1,] <- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/alpha.csv")))
  
  # beta.smpl[1,]  <-  rgamma(I, shape =1 , rate = lambda.beta)
  beta.smpl[1,] <- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/beta.csv")))
  
  
  
  # Initial values of proposal distributions
  
  theta.c <-  theta.smpl[1,]
  alpha.c <-  alpha.smpl[1,]
  beta.c  <-  beta.smpl[1,]
  
  
  ##################################################################
  #           main loop
  ##################################################################
  
  # to evaluate the time
  times <- Sys.time()
  # the loop
  
  for (i in 2:mcmc_size ){
    
    
    
    #########################################################################
    # latent traits
    #########################################################################
    
    ######################################
    # compute the numerator
    ######################################
    
    # the prior of each theta is kumara(theta, alpha.theta, beta.theta)
    
    # the proposal values to the  latent trait candidates since kumar
    # shape1 <- tao.theta
    # shape2 (vector) where theta.c is median
    # shape2 <- log(0.5)/log(1-theta.c^(tao.theta))
    
    
    theta.c<-rkumar(N,tao.theta,log(0.5)/log(1-theta.smpl[i-1, ]^(tao.theta)))
    while(theta.c == 0 || theta.c ==1){
      theta.c<-rkumar(N,tao.theta,log(0.5)/log(1-theta.smpl[i-1, ]^(tao.theta)))
    }
    
    
    aux <- matrix(theta.c,N,I,byrow=FALSE)
    ####alpha.theta.smpl o alpha.theta.c??
    L.post.num <- dkumar(theta.c, alpha.theta.smpl, beta.theta.smpl, log = T)
    + apply(ifelse (x, log(pkumar(aux, matrix(alpha.smpl[i-1,], N, I, byrow = T), 
                                  matrix(beta.smpl[i-1,], N, I, byrow = T))), 
                    log(1 - pkumar(aux, alpha.smpl[i-1,], beta.smpl[i-1,]))),1,sum)
    + dkumar(theta.smpl[i-1,], tao.theta, log(0.5)/log(1-theta.c^(tao.theta)), log = T)
    
     # print("L.post.num")
     # print(L.post.num)
    
    ######################################
    # compute the denominator
    ######################################
    
    aux <- matrix(theta.smpl[i-1,],N,I,byrow=FALSE)
    # print("aux")
    # print(aux)
    
    L.post.den <- dkumar(theta.smpl[i-1,], alpha.theta.smpl, beta.theta.smpl, log = T)
    + apply(ifelse (x, log(pkumar(aux, matrix(alpha.smpl[i-1,], N, I, byrow = T), 
                                  matrix(beta.smpl[i-1,], N, I, byrow = T))), 
                    log(1 - pkumar(aux, alpha.smpl[i-1,], beta.smpl[i-1,]))),1,sum)
    + dkumar(theta.c, tao.theta, log(0.5)/log(1-theta.smpl[i-1,]^(tao.theta)), log = T)
    
    # print("L.post.den")
    # print(L.post.den)
    
    
    # alpha
    # maybe, an extra control to avoid division by cero could be neceesary here
    alpha <-  exp(L.post.num - L.post.den)
    # print("alpha")
    # print(alpha)
    
    alpha <- ifelse(alpha>1,1,alpha)
    # test to accept or reject
    
    accept <- ifelse(alpha>runif(N),TRUE,FALSE)
    # print("accept")
    # print(accept)
    
    accept.rate.t <- accept.rate.t + accept
    theta.smpl[i,] <-ifelse(accept, theta.c, theta.smpl[i-1,])
    # print("theta.smpl[i,]")
    # print(theta.smpl[i,])
    
    #########################################################################
    # item parameters
    ###################################
    ######################################
    # proposal of alpha is trunnormal(alpha.smpl[i-1, ], tao.alpha)
    # proposal of beta is truncnorm(beta.smpl[i-1,], tao.beta)
    #
    # priors
    # prior of alpha is gamma with shape = mu^2/ var and rate = mu / var
    # mu = 1 and var = var.alpha
    ######################################
    # compute the numerator
    ######################################
    
    
    # 0. proposal values
    alpha.c <- rtruncnorm(I, a = 1, b = Inf, alpha.smpl[i-1,], tao.alpha)
    
    beta.c <- rtruncnorm(I, a = 1, b = Inf, beta.smpl[i-1,], tao.beta)
    
    # 
    aux <- matrix(theta.smpl[i,],N,I,byrow=FALSE)
    
    L.post.num <- dgamma(alpha.c,  shape = 1/var.alpha, rate = 1/var.alpha, log = T)
    + dgamma(beta.c, shape = 1/var.beta, rate = 1/var.beta, log = T)
    + apply (ifelse(x, 
                    log(pkumar(aux, matrix(alpha.c, N, I, byrow = T), matrix(beta.c, N, I, byrow = T))),
                    log(1- pkumar(aux, matrix(alpha.c, N, I, byrow = T), matrix(beta.c, N, I, byrow = T)))), 
             2, sum)
    + log(dtruncnorm(alpha.smpl[i-1,], 1, Inf, alpha.c, tao.alpha))
    + log(dtruncnorm(beta.smpl[i-1,], 1, Inf, beta.c, tao.beta))
    
    
    
    
    ######################################
    # compute the denominator
    ######################################
    # 
    # aux = matrix(theta.smpl[i,],N,I,byrow=FALSE)
    
    L.post.den <- dgamma(alpha.smpl[i-1, ],  shape = 1/var.alpha, rate = 1/var.alpha, log = T)
    + dgamma(alpha.smpl[i-1,], shape = 1/var.beta, rate = 1/var.beta, log = T)
    + apply (ifelse(x, 
                    log(pkumar(aux, matrix(alpha.smpl[i-1, ], N, I, byrow = T), matrix(beta.smpl[i-1, ], N, I, byrow = T))),
                    log(1-pkumar(aux, matrix(alpha.smpl[i-1, ], N, I, byrow = T), matrix(beta.smpl[i-1, ], N, I, byrow = T)))),
             2, sum)
    + log(dtruncnorm(alpha.c, 1, Inf, alpha.smpl[i-1, ], tao.alpha))
    + log(dtruncnorm(beta.c, 1, Inf, beta.smpl[i-1,], tao.beta))
    
    
  
    # alpha
    # control division by zero P close zero
    alpha <- exp(L.post.num - L.post.den)
    alpha <-  ifelse(alpha>1,1,alpha)
    
   
    
    # test to accept or reject
    accept <- ifelse(alpha>runif(I),TRUE,FALSE)
    accept.rate.i <- accept.rate.i + accept
    
    
    
    alpha.smpl[i,]<-ifelse(accept, alpha.c, alpha.smpl[i-1,])
    beta.smpl[i,]<- ifelse(accept, beta.c, beta.smpl[i-1,])
    
    # show the iteration
    if(show.iteration){
      cat( i,"\n " )
    }
  }# end for i
  
  
  ####################################
  # show the time
  #######################################
  times = Sys.time() - times
  
  cat( "\n mcmc cicle time: ",times,"\n")
  
  ####################################
  # return the data
  #######################################
  
  accept.rate.t <- accept.rate.t/(mcmc_size)
  accept.rate.i  <- accept.rate.i/(mcmc_size*I)
    l <-  list(theta.smpl=theta.smpl,
             alpha.smpl = alpha.smpl, 
             beta.smpl = beta.smpl, 
             accept.rate.t= accept.rate.t, 
             accept.rate.i = accept.rate.i)
  l
  
}# end irt.Metropolis

# sink()
salida <- KKModel(x, 1000, show.iteration = F, 
                  tao.theta = 10, 
                  var.alpha = 1000, var.beta = 1000,
                  tao.alpha = 50, tao.beta = 50
                )
salida


salida$accept.rate.i


# Recuperación de parámetros

N<-dim(salida$theta)[2]
theta <- rep(0, N)
for ( i in 1:N){
  theta[i]<- mean(salida$theta.smpl[,i])
}
ThetaPobl<- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/theta.csv")))
dif<- max(theta - ThetaPobl)
dif


I<-dim(salida$alpha.smpl)[2]
alpha <- rep(0, I)
for ( i in 1:I){
  alpha[i]<- mean(salida$alpha.smpl[,i])
}
alphaPobl<- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/alpha.csv")))
difAlpha<- max(alpha - alphaPobl)
difAlpha



beta <- rep(0, I)
for ( i in 1:I){
  beta[i]<- mean(salida$beta.smpl[,i])
}
betaPobl<- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/beta.csv")))
difBeta<- max(beta - betaPobl)
difBeta
