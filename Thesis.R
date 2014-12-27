

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


rm(list = ls())

x<-read.table(file="/home/mirt/Documentos/Thesis/R/Simulaciones/1_kstest.txt" ,header=T,sep="")
x<-as.matrix(x)

sink("/home/mirt/Documentos/Thesis/R/Simulaciones/Salida.txt")

KKModel = function (x, mcmc_size=100, show.iteration=TRUE,
                    lambda.alpha.theta = 0.2,
                    lambda.beta.theta = 0.2,
                    
                    tao.alpha.theta = 0.25, 
                    tao.beta.theta = 0.25,
                    
                    tao.theta = 0.2, 
                    tao.alpha = 0.25, 
                    tao.beta=0.2, 
                    
                    lambda.alpha = 0.1, 
                    lambda.beta = 0.1){

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
alpha.theta.smpl <- matrix(0 , nrow = mcmc_size, ncol =1)
beta.theta.smpl <- matrix(0 , nrow = mcmc_size, ncol =1)
alpha.smpl <- matrix(1,nrow = mcmc_size, ncol = I)
 beta.smpl <- matrix(0,nrow = mcmc_size, ncol = I)
theta.smpl <- matrix(0,nrow = mcmc_size, ncol = N)

# Other global values
# aceptation rate fot hyperparameter
accept.rate.h <- 0
# aceptation rate for thetas
accept.rate.t <- 0
# aceptation rate for items
accept.rate.i <- 0


##################################################################
#           initial values
##################################################################
# First values in the chains

# alpha.theta.smpl[1] <-  rgamma(1, shape = 1, rate = lambda.alpha.theta)+1
  alpha.theta.smpl[1] <-  2.315169

# beta.theta.smpl[1]  <-  rgamma(1, shape = 1, rate = lambda.beta.theta)+1
  beta.theta.smpl[1] <- 2.557957

# theta.smpl[1,] <-  rkumar(N, 2, 2.5)
  theta.smpl[1,] <- c(0.5255650, 0.6086310, 0.7127212, 0.1388326, 0.6998396, 0.6875798, 0.8943026, 0.5486967, 0.6072473, 0.5092753, 
                      0.4812270, 0.8025894, 0.7803526, 0.3986964, 0.4242609)
  

# alpha.smpl[1,] <-   rgamma(I, shape =1 , rate = lambda.alpha)
  alpha.smpl[1,] <- c( 5.756020, 7.187109, 6.411305, 3.338797, 8.248473, 7.532518, 1.968525, 7.761372)
  
# beta.smpl[1,]  <-  rgamma(I, shape =1 , rate = lambda.beta)
  beta.smpl[1,] <- c(218.423876, 2181.820969, 550.650100, 1.510153, 9500.053667, 19.213044, 1.917371, 225.858270)
  


# Initial values of proposal distributions
alpha.theta.c <- alpha.theta.smpl[1]
beta.theta.c <- beta.theta.smpl[1]
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
  # hyperparameters 
  #########################################################################

  
  # proposal 
  
  alpha.theta.c <- rtruncnorm(1, 1, Inf, alpha.theta.smpl[i-1], tao.alpha.theta)
  
  beta.theta.c <- rtruncnorm(1, 1, Inf, beta.theta.smpl[i-1], tao.beta.theta)
  
  ######################################
  # compute the numerator
  ######################################
  
  # the prior is exponential with lambda = 0.2
  
  L.post.num =  dexp(alpha.theta.c , lambda.alpha.theta, log = TRUE)
  + dexp(beta.theta.c, lambda.beta.theta, log = TRUE)
  + sum(dkumar(theta.smpl[i-1], alpha.theta.c, beta.theta.c, log = TRUE))
  + log(dtruncnorm(alpha.theta.smpl[i-1], 1, Inf, alpha.theta.c, tao.alpha.theta))
  + log(dtruncnorm(beta.theta.smpl[i-1], 1, Inf, beta.theta.c, tao.beta.theta))
  
  ######################################
  # compute the denominator
  ######################################
  
  L.post.den = dexp(alpha.theta.smpl[i-1], lambda.alpha.theta, log = TRUE)
  + dexp(beta.theta.smpl[i-1], lambda.beta.theta, log = T)
  + sum(dkumar(theta.smpl[i-1], alpha.theta.smpl[i-1], beta.theta.smpl[i-1], log = T))
  + log(dtruncnorm(alpha.theta.c, 1, Inf, alpha.theta.smpl[i-1], tao.alpha.theta))
  + log(dtruncnorm(beta.theta.c, 1, Inf, beta.theta.smpl[i-1], tao.beta.theta))
  
  # 
  alpha <- exp(L.post.num - L.post.den)
  alpha <- ifelse(alpha>1,1,alpha)
  
  # test to accept or reject
  accept <- ifelse(alpha>runif(2),TRUE,FALSE)
  accept
  accept.rate.h <- accept.rate.h + sum(accept)
  alpha.theta.c = alpha.theta.smpl[i] =  ifelse(accept[1], alpha.theta.c, alpha.theta.smpl[i-1])
  beta.theta.c = beta.theta.smpl[i] =  ifelse(accept[2], beta.theta.c, beta.theta.smpl[i-1])
  
  print("alpha.theta.smpl[i]")
  print(alpha.theta.smpl[i])
  print("beta.theta.smpl[i] ")
  print(beta.theta.smpl[i])
  
  #########################################################################
  # latent traits
  #########################################################################
  
  ######################################
  # compute the numerator
  ######################################
  # the prior of each theta is kumara(theta, alpha.theta, beta.theta)
  
  # beta.theta.c (vector) where theta.c is median
  # alpha same all items
    print("theta.smpl[i-1,]")
    print(theta.smpl[i-1,])
    beta.theta.c <- log(0.5)/log(1-theta.smpl[i-1,]^(tao.theta))
    print("beta.theta.c")
    print(beta.theta.c)
    
  # 0. proposal values to the  latent trait candidates: K(alpha = tao,theta, beta = beta.theta.c)
   
  theta.c<-rkumar(length(beta.theta.c),tao.theta,beta.theta.c)
  while(theta.c == 0 || theta.c ==1){
    theta.c<-rkumar(length(beta.theta.c),tao.theta,beta.theta.c)
  }
  
  print("theta.c")
  print(theta.c)
  #
  # 1. 
  

  aux <- matrix(theta.c,N,I,byrow=FALSE)
  
  L.post.num <- dkumar(theta.c, alpha.theta.smpl[i], beta.theta.smpl[i], log = T)
  + apply(ifelse (x, log(pkumar(aux, alpha.theta.smpl[i], beta.theta.smpl[i])), 
                  log(1 - pkumar(aux, alpha.theta.smpl[i], beta.theta.smpl[i]))),1,sum)
  + dkumar(theta.c, tao.theta, beta.theta.c, log = T)

  
    
  ######################################
  # compute the denominator
  ######################################
  # 1. 
  aux = matrix(theta.smpl[i-1,],N,I,byrow=FALSE)
  # print("aux")
  # print(aux)
  
  L.post.den = dkumar(theta.smpl[i-1,], alpha.theta.smpl[i], beta.theta.smpl[i], log = T)
  + apply(ifelse (x, log(pkumar(aux, alpha.theta.smpl[i], beta.theta.smpl[i])), 
                  log(1 - pkumar(aux, alpha.theta.smpl[i], beta.theta.smpl[i]))),1,sum)
  + dkumar(theta.smpl[i-1,], tao.theta, beta.theta.c, log = T)
  
  # print("L.post.den")
  # print(L.post.den)
  
  
  # alpha
  # maybe, an extra control to avoid division by cero could be neceesary here
  alpha = exp(L.post.num - L.post.den)
  alpha = ifelse(alpha>1,1,alpha)
  # test to accept or reject
  accept = ifelse(alpha>runif(N),TRUE,FALSE)
  accept.rate.t = accept.rate.t + sum(accept)
  theta.c = theta.smpl[i,] = ifelse(accept, theta.c, theta.smpl[i-1,])
  print("theta.smpl[i,]")
  print(theta.smpl[i,])
  
  #########################################################################
  # item parameters
  ###################################
  ######################################
  # proposal of alpha is normal(alpha.smpl[i-1, ], tao.alpha)
  # proposal of beta is gamma(beta.smpl[i-1,], tao.beta)
  #
  # priors
  # prior of alpha is exponencial(lambda.alpha = 0.1)
  # prior of beta  is exponencial(lambda.beta = 0.1)
  ######################################
  # compute the numerator
  ######################################
     
    
  # 0. proposal values
  alpha.c = rtruncnorm(I, a = 0, b = Inf, alpha.smpl[i-1,], tao.alpha)
  
  beta.c = rtruncnorm(I, a = 0, b = Inf, beta.smpl[i-1,], tao.beta)

  
    
  # 
  aux = matrix(theta.smpl[i-1,],N,I,byrow=FALSE)
  alpha.matrix = matrix(alpha.c, N, I, byrow = TRUE)
  P = aux^alpha.matrix
  beta.matrix = matrix(beta.c, N, I, byrow = TRUE)
  P = 1- (1- P)^beta.matrix
  
  # 2. log of probability of correct response(KK model)
  L.post.num = matrix(apply(ifelse(x,log(P),log(1-P)),2,sum),I,1)
  
  # 3. log posterior  of each item parameter
  L.pos.num = L.post.num 
  - matrix(lambda.alpha*alpha.c,I,1)- matrix(lambda.beta*beta.c,I,1)
  
 
   ######################################
   # compute the denominator
   ######################################
   # 
   aux = matrix(theta.smpl[i-1,],N,I,byrow=FALSE)
   alpha.matrix = matrix(alpha.smpl[i-1,], N, I, byrow = TRUE)
   P = aux^alpha.matrix
   beta.matrix = matrix(beta.smpl[i-1,], N, I, byrow = TRUE)
   P = 1- (1- P)^beta.matrix

   # 2. log of probability of correct response(KK model)
   L.post.den = matrix(apply(ifelse(x,log(P),log(1-P)),2,sum),I,1)

  # 3. log posterior  of each item parameter
  L.pos.den = L.post.den - 
    matrix(lambda.alpha*alpha.smpl[i-1,],I,1)- matrix(lambda.beta*beta.smpl[i-1,],I,1)

 
  # alpha
  # control division by zero P close zero
  alpha = exp(L.post.num - L.post.den)
  alpha = ifelse(alpha>1,1,alpha)
  # test to accept or reject
  accept = ifelse(alpha>runif(I),TRUE,FALSE)
  accept.rate.i = accept.rate.i + sum(accept)    

  alpha.c = alpha.smpl[i,] =  ifelse(accept, alpha.c, alpha.smpl[i-1,])
  beta.c = beta.smpl[i,] =  ifelse(accept, beta.c, beta.smpl[i-1,])
  
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
  accept.rate.h =accept.rate.h/(mcmc_size)
  accept.rate.t =accept.rate.t/(mcmc_size*N)
  accept.rate.i =accept.rate.i/(mcmc_size*I)
  l = list(theta.smpl=theta.smpl, alpha.smpl = alpha.smpl,   beta.smpl = beta.smpl, accept.rate.h = accept.rate.h, accept.rate.t= accept.rate.t, accept.rate.i=accept.rate.i)
  l

}# end irt.Metropolis


salida = KKModel(x, 1000)
salida

sink()



