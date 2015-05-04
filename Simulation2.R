require(tmvtnorm)
require(pracma)
require(mvtnorm)

######################
# Defining functions #
######################

# Epanechnikov kernel
K = function(t){
  res = 0.75*(1-t^2)*(abs(t)<=1)
  return(res)
}

#####################
# Solving eqn 3 & 6 #
#####################

# evaluating expressions independent of mu and sigma2
intialize = function(data, t_index, theta0, tau, N=1000){
  n = dim(data)[1]
  m = dim(data)[2]
  p = length(theta0)
  
  # initialization values for mu and sigma2
  x_spline = apply(data, 1, splinefun, x=t_index) # length n
  mu_0 = function(t){ # mu_0 to be evaluated at t (vector)
    x_t = matrix(nrow=n, ncol=length(t))
    for (i in 1:n){
      x_t[i,] = x_spline[[i]](t)
    }
    return(colMeans(x_t))
  }
  sigma2_0 = sum((data-mean(data))^2)/(n*m)
  
  # return N random sample of theta
  theta_samples = matrix(NA, N, p)
  i = 0
  while (i<N){
    samp = rtmvnorm(n=1, mean=theta0, sigma=diag(tau^2,p,p), # N-by-p
                    lower=rep(min(t_index),p), upper=rep(max(t_index),p), 
                    algorithm='rejection')
    samp_order = sum((order(samp)==1:p))
    if (samp_order==p){
      i = i+1
      theta_samples[i,] = samp
    }
  }
  
  # g(t_ij, theta_l), independent of i for standardized sampling time
  g_mat = t(apply(theta_samples, 1, pchip, yi=theta0, x=t_index)) # N-by-m
  # s_ij, independent of i for standardized sampling time
  s_vec = apply(g_mat, 2, sd)
  # lambda for kernel
  lambda = (((243*3/5)/(35*N/25))^0.2)*mean(s_vec)
  return(list(mu_0=mu_0, sigma2_0=sigma2_0, 
              theta_samples=theta_samples, g_mat=g_mat, lambda=lambda))
}

# evaluating mu and sigma2 till convergence
eval = function(data, t_index, init_res, tol=10^-4, max_itr = 1000){
  mu_0 = init_res$mu_0
  sigma2_0 = init_res$sigma2_0
  theta_samples = init_res$theta_samples
  g_mat = init_res$g_mat
  lambda = init_res$g_mat
  
  tol_mu = rep(NA, max_itr)
  tol_sigma2 = rep(NA, max_itr)
  
  n = dim(data)[1]
  m = dim(data)[2]
  N = dim(theta_samples)[1]
  p = length(theta0)
  
  itr = 0
  mu_old = splinefun(x=t_index, y=t_index) # distance evaluated as L2 norm at t_index
  sigma2_old = 100
  if (itr==0){
    mu_new = mu_0
    sigma2_new = sigma2_0
  }
  
  while (itr<max_itr & (abs(sigma2_new-sigma2_old)>tol |
                          sqrt(sum((mu_new(t_index)-mu_old(t_index))^2))>tol)){
    itr = itr+1
    tol_mu[itr] = sqrt(sum((mu_new(t_index)-mu_old(t_index))^2))
    tol_sigma2[itr] = abs(sigma2_new-sigma2_old)
    print(itr)
    mu_old = mu_new
    sigma2_old = sigma2_new
    # mu(g(t_i,theta_l)), independent of i for standardized sampling time
    mu_g = t(apply(g_mat, 1, mu_old)) # N-by-m
    # f(x_i|theta_l)
    fx_theta = matrix(NA, N, n) # N-by-n
    for (i in 1:n){
      fx_theta[,i] = apply(mu_g, 1, dmvnorm, x=data[i,], sigma=diag(sigma2_old,m,m))
    }
    # est f(x_i) 
    fx = colMeans(fx_theta) # length n
    # est w_ij(s|x_i)
    w = list()
    for (i in 1:n){
      w[[i]] = list()
      for (j in 1:m){
        fg_x = function(s){
          summand = 0
          for (l in 1:N){
            summand = summand+K((g_mat[l,j]-s)/lambda)*fx_theta[l,i]/fx[i]
          }
          return(summand/(N*lambda))
        }
        w[[i]][[j]] = fg_x
      }
    }
    
    # calculate mu_new
    mu_new = function(s){
      numerator = 0
      denominator = 0
      for (i in 1:n){
        for (j in 1:m){
          numerator = numerator+data[i,j]*w[[i]][[j]](s)
          denominator = denominator + w[[i]][[j]](s)
        }
      }
      return(numerator/denominator)
    } 
    
    # calculate sigma2_new
    sum_inner = matrix(NA, N, n)
    for (l in 1:N){
      for (i in 1:n){
        sum_inner[l,i] = sum((data[i,]-mu_g[l,])^2)*fx_theta[l,i]/fx[i]
      } 
    }
    sigma2_new = sum(colMeans(sum_inner))/(n*m)
  }
  
  # log-likelihood
  loglike = sum(log(fx))
  
  return(list(mu=mu_new, mu_vec = mu_new(t_index), sigma2=sigma2_new, loglike=loglike,
              tol_mu = tol_mu, tol_sigma2 = tol_sigma2))
}

##################
# Implementation #
##################

# data, n-by-m_i (m_i = m in our case)
# setwd('/home/wenyu/Desktop/Cornell/Class/ORIE 6070/Project/R Code')
subj0 = read.csv(file="subj0.csv",head=FALSE,sep=",") # 53-by-128
data = as.matrix(subj0)
t_index = (1:128)*2
theta0 = c(40,120,160,240)
tau = 3

init_res = intialize(data, t_index, theta0, tau, N=1000)
eval_res = eval(data, t_index, init_res, tol=0.01, max_itr = 1000)

prefix = 'subj0_'
write.csv(eval_res$mu_vec, file = paste(prefix, "mu_vec.csv", sep=''))
write.csv(eval_res$sigma2, file = paste(prefix, "sigma2.csv", sep=''))
write.csv(eval_res$loglike, file = paste(prefix, "loglike.csv", sep=''))
write.csv(eval_res$tol_mu, file = paste(prefix, "tol_mu.csv", sep=''))
write.csv(eval_res$tol_sigma2, file = paste(prefix, "tol_sigma2.csv", sep=''))

sink(paste(prefix, ".txt", sep=''), append=FALSE, split=FALSE)
