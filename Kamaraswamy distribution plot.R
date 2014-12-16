
alpha=seq(from=0,to=10,by=0.01)
w = 0.9
fk = function(alpha,theta){
  
  beta = log(0.5)/log(1-w^(alpha))
  alpha *  beta * theta^(alpha-1) * (1-theta^alpha)^(beta-1)
}

alpha = 5
beta = 10

theta = 0.9
y = fk(alpha=20,theta=theta)
y
y = fk(alpha=alpha, w =w)
y


plot(alpha,y,type="l",col="blue",main="Kuramaswamy density")
lines(x,y,col="red")

Fk = function(x,alpha,beta){
  1-(1-x^alpha)^beta
}

Y = Fk(alpha=alpha,w =w)

plot(alpha,Y,type="l",col="blue",main="Kuramaswamy distriution")


################################################################

x = seq(from = 0, to = 1, by = 0.01)
y = dkumar(x, 5, 2)
plot(x, y, type = "l")



x = seq(from = -3, to = 3, by = 0.01)
y = dnorm(x, 2, 1)
plot(x, y, type = "l")


############################################################

# alpha vs difficulty

dfk <- function (alpha, w ,theta){
(w^alpha)*log (0.5)*log (w)/(1- w^alpha)*
    ((1-theta^alpha)^(log (0.5)/log (1- w^alpha)-1))*alpha * theta^(alpha -1)+
    log (0.5)/log(1-w^alpha)*((theta^(alpha -1 )+ alpha*theta^(alpha -1))* 
                              ((1-theta^alpha)^(log (0.5)/log (1- w^alpha)-1))+
                              ((1-theta^alpha)^(log (0.5)/log (1- w^alpha)-1))*
                              ((w^alpha)*log (0.5)*log (w)/(1- w^alpha)-
                                 (log (0.5)/log(1-w^alpha)-1)*(theta^alpha)*log (theta)/ (1- theta^alpha))*
                                 alpha*(theta^(alpha -1)) )
                            
}

# Proof with constant theta 
w = 0.7

theta = 0.1
alpha <- seq(from = 0.01, to = 10, by = 0.1)
median<-c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
#windows()
plot(alpha,dfk(alpha,.1,theta),type="l", ylim = c(0, 30))
for(i in 1:length(median)){
  points(alpha,dfk(alpha,median[i], theta),type="l",col=rainbow(length(median))[i])  
}
legend(6,30,legend=median,col=rainbow(length(median)),lty=1)

# Proof with constant w


w = 0.8
alpha <- seq(from = 1, to = 10, by = 0.1)
theta<-c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
#windows()
plot(alpha,dfk(alpha,w,.1),type="l", ylim = c(0, 30))
for(i in 1:length(theta)){
  points(alpha,dfk(alpha, w, theta[i]),type="l",col=rainbow(length(theta))[i])  
}
legend(6,30,legend=median,col=rainbow(length(theta)),lty=1)


