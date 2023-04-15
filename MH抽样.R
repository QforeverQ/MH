N = 3000  #总抽样次数
n = 100 #样本个数
set.seed(1234)
x = rnorm(n, 200, sqrt(2))  #从给定分布中抽取100个样本

#初始值
mu0 = 0
phi0 = 5

#输出结果
para = matrix(0,N,2)
para[1, ]=c(mu0, phi0)

#Metropolis-Hastings算法
metropolis<-function(x,mu,phi,n=100){
  mustar = rnorm(1, mu, 1)
  rat = exp((-phi/2)*(n*(mustar^2-mu^2)-2*sum(x)*(mustar-mu))) * 
        (dnorm(mu,mustar,1)/dnorm(mustar,mu,1))
  if (runif(1) > rat)  mustar = mu
  phistar = rnorm(1, phi, 0.1)
  rat = ((phistar/phi)^(n/2-1)*exp((-(phistar-phi)/2)*sum((x-mustar)^2))) *
        (dnorm(phi,phistar,0.1)/dnorm(phistar,phi,0.1))
  if (runif(1) > rat)  phistar = phi
  return (c(mustar,phistar))
}
for (i in 2:N){
  para[i, ] = metropolis(x,para[i-1,1],para[i-1,2])
}

#预烧次数
m = 500 

#图形展示
plot(para[,1], type='l', xlab = NULL, ylab = "All tracks of mu")
plot(para[,1], type='l', xlim = c(580,2920), ylim = c(199,200.5), 
     ylab = "Burn-in tracks of mu", xaxt = 'n')
axis(1,c(500,1000,1500,2000,2500,3000))
plot(para[,2], type='l', xlab = NULL, ylab = "All tracks of phi")
plot(para[,2], type='l', xlim = c(580,2920), ylim = c(0.2,0.8), 
     ylab = "Burn-in tracks of phi", xaxt = 'n')
axis(1,c(500,1000,1500,2000,2500,3000))
hist(para[m:N,1], prob=TRUE, xlab = "mu")
hist(para[m:N,2], prob=TRUE, xlab = "phi")
acceptmu = length(unique(para[,1]))/N  #mu的接受率=0.255 
acceptphi = length(unique(para[,2]))/N #phi的接受率=0.541
mean(para[m:N,1])  #mu的采样均值=199.782
mean(para[m:N,2])  #phi的采样均值=0.498


#Metropolis-Hastings-Gibbs
for (i in 2:N){
  para[i,1] = rnorm(1, sum(x)/n, 1/sqrt(n*para[i-1,2]))
  para[i,2] = rgamma(1, n/2, sum((x-para[i,1])^2)/2)
}

#预烧次数
m = 10 

#图形展示
plot(para[,1], type='l', xlab = NULL, ylab = "All tracks of mu")
plot(para[2:3000,1], type='l', xlim = c(80,2920), ylim = c(199,200.5), 
     ylab = "Burn-in tracks of mu", xaxt = 'n')
axis(1,c(0,500,1000,1500,2000,2500,3000))
plot(para[,2], type='l', xlab = NULL, ylab = "All tracks of phi")
plot(para[2:3000,2], type='l', xlim = c(80,2920), ylim = c(0.2,0.8), 
     ylab = "Burn-in tracks of phi", xaxt = 'n')
axis(1,c(0,500,1000,1500,2000,2500,3000))
hist(para[2:N,1], prob=TRUE, xlab = "mu")
hist(para[2:N,2], prob=TRUE, xlab = "phi")
mean(para[2:N,1])  #mu的采样均值=199.777
mean(para[2:N,2])  #phi的采样均值=0.496
