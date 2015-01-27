###################################
# Universidad Nacional de Colombia
# Yuli Chacon
# Director A. Montenegro
# Algoritmo MH
# Modelo KK
# Estimación de los trazos, parámetros de ítems fijos.
###################################

# utilities libraries


library("VGAM")
library("truncnorm")
library(car)
library(coda)


rm(list = ls())

x<-read.table(file="/home/mirt/Documentos/Thesis/R/Simulaciones/1_kstest.txt" ,header=T,sep="")
x<-as.matrix(x)

# sink("/home/mirt/Documentos/Thesis/R/Simulaciones/SalidaOnlyTraits.txt")

KKModel = function (x, mcmc_size=2, show.iteration=TRUE,
                    var.alpha.theta = 100,
                    var.beta.theta = 100,
                    
                    tao.alpha.theta = 0.5, 
                    tao.beta.theta = 0.5,
                    
                    tao.theta = 1 ){
  
  # x: test data
  # mcmc_size: number of samples
  # file.name: file name to save the results, recomended)
  # show.iteration=TRUE: show each iteration
  # lambda.alpha.theta: prior distribution parameter hyperparemeter alpha.theta
  # lambda.beta.theta: prior distribution parameter hyperparemeter beta.theta 
  # alpha.theta: prior distribution parameter theta
  # beta.theta: prior distribution parameter theta
  # tao.theta: d  tunning parameter for latent traits (proposal)
  
  
  
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
  alpha.theta.smpl <- matrix(0 , nrow = mcmc_size, ncol =1)
  beta.theta.smpl <- matrix(0 , nrow = mcmc_size, ncol =1)
  theta.smpl <- matrix(0,nrow = mcmc_size, ncol = N)
  
  # Other global values
  # aceptation rate fot hyperparameter
  accept.rate.h <- c(0,0)
  # aceptation rate for thetas
  accept.rate.t <- rep(0, N)

  
  ##################################################################
  #           initial values
  ##################################################################
  # First values in the chains
  
  # alpha.theta.smpl[1] <-  rgamma(1, shape = 1, rate = lambda.alpha.theta)+1
  alpha.theta.smpl[1] <-   4.938939
  
  # beta.theta.smpl[1]  <-  rgamma(1, shape = 1, rate = lambda.beta.theta)+1
  beta.theta.smpl[1] <- 6.696828
  
  
  # theta.smpl[1,] <-  rkumar(N, 2, 2.5)
  theta.smpl[1,] <- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/theta.csv")))

  
  # alpha.smpl <-   rgamma(I, shape =1 , rate = lambda.alpha)
  alpha.smpl<- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/alpha.csv")))
  
  # beta.smpl  <-  rgamma(I, shape =1 , rate = lambda.beta)
  beta.smpl <- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/beta.csv")))
  
  
  # Initial values of proposal distributions
  
  alpha.theta.c <- alpha.theta.smpl[1]
  beta.theta.c <- beta.theta.smpl[1]
  theta.c <-  theta.smpl[1,]
  
  
  ##################################################################
  #           main loop
  ##################################################################
  
  # to evaluate the time
  times <- Sys.time()
  # the loop
  
  for (i in 2:mcmc_size ){
    
    #########################################################################
    # hyperparameters 
    #########################################################################
    
    
    # proposal 
    
    alpha.theta.c <- rtruncnorm(1, 1, Inf, alpha.theta.smpl[i-1], tao.alpha.theta)
    
    beta.theta.c <- rtruncnorm(1, 1, Inf, beta.theta.smpl[i-1], tao.beta.theta)
    
    ######################################
    # compute the numerator
    ######################################
    
    # the prior is gamma with alpha = mu^2 / var, beta = mu^2 / var
    # mu = 1, var = var.alpha.theta
    
    L.post.num.alpha <- dgamma(alpha.theta.c, shape = 1/var.alpha.theta, rate = 1/var.alpha.theta, log = T)
    + sum(dkumar(theta.smpl[i-1], alpha.theta.c, beta.theta.c, log = TRUE))
    + log(dtruncnorm(alpha.theta.smpl[i-1], 1, Inf, alpha.theta.c, tao.alpha.theta))
    
    
    L.post.num.beta <- dgamma(beta.theta.c,  shape = 1/var.beta.theta, rate = 1/var.beta.theta, log = T)
    + sum(dkumar(theta.smpl[i-1], alpha.theta.c, beta.theta.c, log = TRUE))
    + log(dtruncnorm(beta.theta.smpl[i-1], 1, Inf, beta.theta.c, tao.beta.theta))
    
    ######################################
    # compute the denominator
    ######################################
    
    L.post.den.alpha <- dgamma(alpha.theta.smpl[i-1],  shape = 1/var.alpha.theta, rate = 1/var.alpha.theta, log = T)
    + sum(dkumar(theta.smpl[i-1], alpha.theta.smpl[i-1], beta.theta.smpl[i-1], log = T))
    + log(dtruncnorm(alpha.theta.c, 1, Inf, alpha.theta.smpl[i-1], tao.alpha.theta))
    
    L.post.den.beta <- dgamma(beta.theta.smpl[i-1],  shape = 1/var.beta.theta, rate = 1/var.beta.theta, log = T)
    + sum(dkumar(theta.smpl[i-1], alpha.theta.smpl[i-1], beta.theta.smpl[i-1], log = T))
    + log(dtruncnorm(beta.theta.c, 1, Inf, beta.theta.smpl[i-1], tao.beta.theta))
    
    alpha <- exp(L.post.num.alpha - L.post.den.alpha)
    alpha <- ifelse(alpha>1,1,alpha)
    
    beta <- exp(L.post.num.beta - L.post.den.beta)
    beta <- ifelse(beta>1,1, beta)
    
    # test to accept or reject
    accept.alpha <- ifelse(alpha>runif(1),TRUE,FALSE)
    accept.rate.h[1] <- accept.rate.h[1] + sum(accept.alpha)
    
    accept.beta <- ifelse(beta>runif(1),TRUE,FALSE)
    accept.rate.h[2] <- accept.rate.h[2] + sum(accept.beta)
    
    alpha.theta.smpl[i] <-ifelse(accept.alpha, alpha.theta.c, alpha.theta.smpl[i-1])
    beta.theta.smpl[i] <-ifelse(accept.beta, beta.theta.c, beta.theta.smpl[i-1])
    
    # print("alpha.theta.smpl[i]")
    # print(alpha.theta.smpl[i])
    # print("beta.theta.smpl[i] ")
    # print(beta.theta.smpl[i])
    
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
    L.post.num <- dkumar(theta.c, alpha.theta.smpl[i], beta.theta.smpl[i], log = T)
    + apply(ifelse (x, log(pkumar(aux, matrix(alpha.smpl, N, I, byrow = T), 
                                  matrix(beta.smpl, N, I, byrow = T))), 
                    log(1 - pkumar(aux, alpha.smpl, beta.smpl))),1,sum)
    + dkumar(theta.smpl[i-1,], tao.theta, log(0.5)/log(1-theta.c^(tao.theta)), log = T)
    
    # print("L.post.num")
    # print(L.post.num)
    
    ######################################
    # compute the denominator
    ######################################
    
    aux <- matrix(theta.smpl[i-1,],N,I,byrow=FALSE)
    # print("aux")
    # print(aux)
    
    L.post.den <- dkumar(theta.smpl[i-1,], alpha.theta.smpl[i], beta.theta.smpl[i], log = T)
    + apply(ifelse (x, log(pkumar(aux, matrix(alpha.smpl, N, I, byrow = T), 
                                  matrix(beta.smpl, N, I, byrow = T))), 
                    log(1 - pkumar(aux, alpha.smpl, beta.smpl))),1,sum)
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
  accept.rate.h <- accept.rate.h/(mcmc_size)
  accept.rate.t <- accept.rate.t/(mcmc_size)
  
  l <-  list(alpha.theta.smpl=alpha.theta.smpl, 
             beta.theta.smpl=beta.theta.smpl, 
             theta.smpl=theta.smpl,
             accept.rate.h = accept.rate.h, 
             accept.rate.t= accept.rate.t
             )
  l
  
}# end irt.Metropolis

# sink()
salida <- KKModel(x, 5, show.iteration = F, 
                  var.alpha.theta = 1000, var.beta.theta = 1000, 
                  tao.alpha.theta = 50,tao.beta.theta = 50, 
                  tao.theta = 10
                  )
salida


salida$accept.rate.i

salida$accept.rate.h

