

###################################
# Universidad Nacional de Colombia
# Yuli Chacon
# Director A. Montenegro
# Algoritmo MH
# Modelo KK
# Estimación de los parámetros de los ítems, habilidades fijas.
###################################

# utilities libraries


library("VGAM")
library("truncnorm")
library(car)
library(coda)


rm(list = ls())

x<-read.table(file="/home/mirt/Documentos/Thesis/R/Simulaciones/1_kstest.txt" ,header=T,sep="")
x<-as.matrix(x)

sink("/home/mirt/Documentos/Thesis/R/Simulaciones/SalidaOnlyItems.txt")

KKModel = function (x, mcmc_size=10, show.iteration=TRUE,
                                                      
                    tao.alpha = 0.25, 
                    tao.beta=0.25, 
                    
                    var.alpha = 100, 
                    var.beta = 100
                    
                    fixed = fixed,
                    alpha.fixed = alpha.fixed,
                    beta.fixed = beta.fixed){
  
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
  
  
  # Other global values

  # aceptation rate for items
  accept.rate.alpha <- rep(0, I)
  accept.rate.beta <- rep(0, I)
  
  ##################################################################
  #           initial values
  ##################################################################
  # First values in the chains
  
  # alpha.theta.smpl <-  rgamma(1, shape = 1, rate = lambda.alpha.theta)+1
  alpha.theta.smpl<-  2.315169
  
  # beta.theta.smpl  <-  rgamma(1, shape = 1, rate = lambda.beta.theta)+1
  beta.theta.smpl <- 2.557957
  
  # theta.smpl <-  rkumar(N, 2, 2.5)
  theta.smpl<- c(0.5255650, 0.6086310, 0.7127212, 0.1388326,0.6998396, 
                      0.6875798,0.8943026, 0.5486967, 0.6072473, 0.5092753, 
                      0.4812270, 0.8025894, 0.7803526, 0.3986964, 0.4242609)
  
  
  # alpha.smpl[1,] <-   rgamma(I, shape =1 , rate = lambda.alpha)
  alpha.smpl[1,] <- c( 5.756020, 7.187109, 6.411305, 3.338797, 8.248473, 7.532518, 1.968525, 7.761372)
  
  # beta.smpl[1,]  <-  rgamma(I, shape =1 , rate = lambda.beta)
  beta.smpl[1,] <- c(218.423876, 2181.820969, 550.650100, 1.510153, 9500.053667, 19.213044, 1.917371, 225.858270)
  
  
  
  # Initial values of proposal distributions
  
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
    
    L.post.num.alpha <- rep(0, I)
    L.post.num.beta <- rep(0, I)
    
    for(l in 1:I){
      if( l != fixed.parameter){
        for( m in 1:N){
        tem.num.alpha <- dgamma(alpha.c[l],  shape = 1/var.alpha, rate = 1/var.alpha, log = T) 
        + ifelse(x[l,m], log(pkumar(theta.smpl[m], alpha.c[l], beta.c[l])), 
                 log(1 - pkumar(theta.smpl[m], alpha.c[l], beta.c[l])) )
        + log(dtruncnorm(alpha.smpl[i-1,l], 1, Inf, alpha.c[l], tao.alpha))
        L.post.num.alpha[l] <- L.post.num.alpha[l] + temp.num.alpha
        
        temp.num.beta <- dgamma (beta.c[l], shape = 1/var.beta, rate = 1/var.beta, log = T )
        + ifelse(x[l,m], log(pkumar(theta.smpl[m], alpha.c[l], beta.c[l])), 
                 log(1 - pkumar(theta.smpl[m], alpha.c[l], beta.c[l])) )
        + log(dtruncnorm(beta.smpl[i-1,l], 1, Inf, beta.c[l], tao.beta)) 
        L.post.num.beta[l] <- L.post.num.beta[l] + temp.num.beta
        
        }
      }
      else{
        L.post.num.alpha[l] <- fixed.alpha[l]
        L.post.num.beta[l] <- fixed.beta[l]
      }
      L.post.num.alpha
      L.post.num.beta
    }
    
    
    ######################################
    # compute the denominator
    ######################################
    
    L.post.den.alpha <- rep(0, I)
    L.post.den.beta <- rep(0, I)
      
    for( l in 1:I){
      if (l != fixed.parameter){
        for (m in 1:N){
          temp.den.alpha <- dgamma(alpha.smpl[i-1,m],  shape = 1/var.alpha, rate = 1/var.alpha, log = T)
          + ifelse(x[l,m], log(pkumar(theta[m], alpha.smpl[i-1, m], beta.smpl[i-1, m])), 
                   log(1- pkumar(theta[m], alpha.smpl[i-1, m], beta.smpl[i-1, m])))
          + log(dtruncnorm(alpha.c[l], 1, Inf, alpha.smpl[i-1,m ], tao.alpha))
          L.post.den.alpha[l] <- L.post.num.alpha[l] + temp.den.alpha
          
          temp.den.beta <- dgamma (beta.smpl[i-1, m], shape = 1/var.beta, rate = 1/var.beta, log = T)
          + ifelse(x[l,m], log(pkumar(theta[m], alpha.smpl[i-1, m], beta.smpl[i-1, m])), 
                   log(1- pkumar(theta[m], alpha.smpl[i-1, m], beta.smpl[i-1, m])))
          + log(dtruncnorm(beta.c[l], 1, Inf, beta.smpl[i-1,m], tao.beta))
          L.post.den.beta[l]<-L.post.den.beta[l] + temp.den.beta
        }
      }
      else {
        L.post.den.alpha<-fixed.alpha[l]
        L.post.den.beta<-fixed.beta[l]
        }
      
      L.post.den.alpha
      L.post.den.beta
    }
  
    # alpha
    # control division by zero P close zero
    
    alpha <-rep(0, I)
    beta <- rep(0, I)
        
    accept.alpha <-rep(0, I)
    accept.beta <- rep(0, I)

    
    for (l in 1:I){
      ifelse( l != fixed){
        
        alpha[l] <- exp(L.post.num.alpha[l] - L.post.den.alpha[l])
        alpha[l] <- ifelse(alpha[l]> 1, 1, alpha[l])
        
        beta[l] <- exp(L.post.num.beta[l] - L.post.den.beta[l])
        beta[l] <- ifelse(beta[l]>1, 1, beta[l])
        
        # test to accept or reject
        # significative change ??
        
        accept.alpha[l] <- ifelse(alpha[l]> runif(1), T, F)
        accept.rate.alpha[l] <- accept.rate.alpha[l] + accept.alpha[l]
       
        accept.beta[l] <- ifelse (beta[l] > runif(1), T, F)
        accept.rate.beta[l] <- accept.rate.beta[l] + accept.beta[l]
        
        alpha.smpl[i,l]<-ifelse(accept.alpha, alpha.c, alpha.smpl[i-1,])
        beta.smpl[i,l]<- ifelse(accept.beta, beta.c, beta.smpl[i-1,])
        
      }
      else{
        alpha.smpl[i,l] <-fixed.alpha
        beta.smpl[i,l] <- fixed.beta
      }
    }
    
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

  accept.rate.alpha <- accept.rate.alpha/(mcmc_size)
  accept.rate.beta <- accept.rate.beta/(mcmc_size)
  l <-  list(alpha.smpl = alpha.smpl, 
             beta.smpl = beta.smpl,
             accept.rate.alpha = accept.rate.alpha,
             accept.rate.beta = accept.rate.beta)
  l
  
}# end irt.Metropolis

sink()
salida <- KKModel(x, 10000, show.iteration = F, 
                  var.alpha = 1000, var.beta = 1000,
                  tao.alpha = 100, tao.beta = 50, 
                  fixed = ,
                  alpha.fixed = ,
                  beta.fixed = )
salida





