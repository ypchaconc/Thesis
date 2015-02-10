

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

# sink("/home/mirt/Documentos/Thesis/R/Simulaciones/SalidaOnlyItems.txt")

  
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
 
  
  ##################################################################
  #              global constants
  ##################################################################
  # N: sample size
  # I: test size
  N <- dim(x)[1]
  I <- dim(x)[2]
  

 
  
  
  ##################################################################
  #           main loop
  ##################################################################
  

    
    
    #########################################################################
    # item parameters
    ###################################
    #
    # priors
    # prior of alpha is gamma with shape = mu^2/ var and rate = mu / var
    # mu = 1 and var = var.alpha
    theta.smpl<- c(t(read.csv(file = "/home/mirt/Documentos/Thesis/R/Simulaciones/theta.csv")))

   
    
    #fx <- function(alpha, beta, var.alpha, var.beta){ 
    fx <- function(param, var.alpha, var.beta,x,theta.smpl){ 
      alpha = param[1]
      beta = param[2]
    dgamma(alpha,  shape = 1/var.alpha, rate = 1/var.alpha) *
      dgamma(beta, shape = 1 /var.beta, rate = 1/var.beta) * 
       sum(ifelse(x, 
                    pkumar(theta.smpl, alpha, beta),
                    1- pkumar(theta.smpl, alpha, beta)))
    }
    
    x0 <- c(runif(1, 1, 50), runif(1, 1, 50) )
    
    

    #Función que ejecuta el slice sampler
    #f es la función a muestrear
    #x0 es el valor inicial dentro del dominio de la función
    #n es el tamaño de la muestra
    #En el ejemplo la función es una normal estándar
    #vamos a suponer que el dominio de la función es el intervalo (-4,4)
    
    slice.sampler = function(f, x0, mcmc_size,var.alpha,var.beta,x, theta.smpl){
      chain = matrix(0,nrow = mcmc_size ,ncol = 2)
      inner.counter = 0
      for (i in 1:mcmc_size){
        y = runif(1,0, f(x0,var.alpha,var.beta,x,theta.smpl))
        found = FALSE
        c = 0
        while (!found){
          x1 <- c(runif(1, 1, 50), runif(1, 1, 50) )
          if(f(x1,var.alpha,var.beta,x,theta.smpl)>y){
            found =TRUE
            chain[i,]<- x1
            x0 <- x1
          }
          c = c+1
        }
        inner.counter <- inner.counter+c
        if(inner.counter %% 1000 == 0){
          print(inner.counter)
        }
      }
      print(paste("Aceptation rate = ",mcmc_size / inner.counter,sep = ""))
      chain
    }
    
    sample = slice.sampler(f = fx,x0 = x0,mcmc_size = 1000,var.alpha = 1000,var.beta = 1000,x,theta.smpl)
  
  
    