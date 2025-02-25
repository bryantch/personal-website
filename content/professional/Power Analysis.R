library(ggplot2)
library(reshape2)
# False Positive Rate Analysis ----
# False Positive Rate - Binomial Distribution ----
N = c(1,2,5,10,25,50,100,150) # vector of interim analyses
alpha = c(.10,.05) # vector of significance levels
p0 = c(.25,.5,.75) # vector of null reliabilities

n = c(1,5,10) # vector of test quantities

seq.binom = function(alpha,p0,n,max){
  Sm = 0
  m = 0
  k = 0
  continue = TRUE
  while(continue) {
    m = m + 1
    k = k + n
    am = qbinom(p = alpha/2, size = k, prob = p0, lower.tail = TRUE)
    bm = qbinom(p = alpha/2, size = k, prob = p0, lower.tail = FALSE)
    
    Xm = rbinom(n = 1, size = n, prob = p0)
    Sm = Sm + sum(Xm)
    if (Sm < am || Sm > bm ||m > max){continue = FALSE}
  }
  return(c(M = m, aM = am, Sum = Sm, bM = bm))
}
check.binom = function(alpha,p0,n,max){
  res = seq.binom(alpha,p0,n,max)
  reject = res[2] > res[3] || res[3] > res[4]
  return(reject)
}
B = 10000
tab.tp = array(data = NA,dim = c(length(N),
                                 length(n),
                                 length(p0),
                                 length(alpha)),dimnames = list(c(N),
                                                                c(n),
                                                                c(p0),
                                                                c(alpha)))

for(l in 1:length(alpha)){
  for(k in 1:length(p0)){
    for(j in 1:length(n)){
      for(m in 1:length(N)){
        tab.tp[m,j,k,l] = sum(replicate(B,check.binom(alpha = alpha[l],
                                                      p0 = p0[k],
                                                      max = N[m],
                                                      n = n[j])))/B
      }
    }
  }
}

tab.tp = tab.tp %>% melt(varnames = c("N","n","p0","alpha"))

labeller = label_bquote(rows = alpha ==. (alpha),
                        cols = p[0] ==. (p0))

tab.tp %>% ggplot() + 
  geom_line(mapping = aes(x = N,y = value,col = as.factor(n))) +
  facet_grid(rows = vars(alpha),cols = vars(p0),labeller = labeller) + ylim(0,1) + 
  theme_bw() + 
  labs(title = "False Positive Rate",
       subtitle = expression(paste("False Positive Rate as a function of p"[0], ", Sample Size (n), # of Interim Analyses (N) and Significance Lvl (",alpha,")")),
       x = expression(paste("# of Interim Analyses (N)")),
       y = expression(paste("P[X" >= "x|p = p"[0],"]")),
       color = "Sample Size (n)")


# False Positive Rate - Normal Distribution ----


# Power Analysis ----
# True Positive Rate - Binomial Distribution ----
N = c(5,10,25,50,100,150) # vector of interim analyses
n = c(1,10) # vector of test quantities

alpha = c(.10,.05) # vector of significance levels
p0 = c(.25,.5,.75) # vector of null reliabilities
pa = c(seq(.10,.30,.10),seq(.40,.60,.05),seq(.70,.90,.10)) # vector of alternative reliabilities

seq.binom = function(alpha,p0,pa,n,max){
  Sm = 0
  m = 0
  k = 0
  continue = TRUE
  while(continue) {
    m = m + 1
    k = k + n
    am = qbinom(p = alpha/2, size = k, prob = p0, lower.tail = TRUE)
    bm = qbinom(p = alpha/2, size = k, prob = p0, lower.tail = FALSE)
    
    Xm = rbinom(n = 1, size = n, prob = pa)
    Sm = Sm + sum(Xm)
    if (Sm < am || Sm > bm ||m > max){continue = FALSE}
  }
  return(c(M = m, aM = am, Sum = Sm, bM = bm))
}

check.binom = function(alpha,p0,pa,n,max){
  res = seq.binom(alpha,p0,pa,n,max)
  reject = res[2] > res[3] || res[3] > res[4]
  return(reject)
}

B = 1000
sum(replicate(B,check.binom(alpha,p0,pa,n,max = 10)))/B

tab.tp = array(data = NA,dim = c(length(pa),
                                 length(N),
                                 length(n),
                                 length(p0),
                                 length(alpha)),dimnames = list(c(pa),
                                                               c(N),
                                                               c(n),
                                                               c(p0),
                                                               c(alpha)))

for(l in 1:length(alpha)){
  for(k in 1:length(p0)){
    for(j in 1:length(n)){
      for(m in 1:length(N)){
        for(i in 1:length(pa)){
          
          tab.tp[i,m,j,k,l] = sum(replicate(B,check.binom(alpha = alpha[l],
                                                        p0 = p0[k],
                                                        max = N[m],
                                                        n = n[j],
                                                        pa = pa[i])))/B
        }
      }
    }
  }
}

library(ggplot2)
library(reshape2)

tab.tp = tab.tp %>% melt(varnames = c("pa","N","n","p0","alpha"))

labeller = label_bquote(rows = alpha ==. (alpha),
                        cols = p[0] ==. (p0))

tab.tp %>% ggplot() + 
  geom_line(mapping = aes(x = pa,y = value,col = as.factor(n),linetype = as.factor(N))) +
  facet_grid(rows = vars(alpha),cols = vars(p0),labeller = labeller) +
  theme_bw() + 
  labs(title = "Power Curves",
       subtitle = expression(paste("Power as a function of p"[0], ", p"[a], ", Sample Size (n) and Significance Lvl (",alpha,")")),
       x = expression(paste("p"[a])),
       y = expression(paste("P[X" >= "x|p = p"[a],"]")),
       linetype = "# of Interim Analyses (N)",
       color = "Sample Size (n)")

# True Positive Rate - Normal Distribution ----
n = c(10,25,50,100,150) # vector of interim analyses
alpha = c(.10,.05,.01) # vector of significance levels
mu_a = c(-2,-1,0,1,2) # vector of alternative means
sigma_a = c(.5,1,2,3,4) # vector of alternative sigmas

seq.norm = function(mu_a,sigma_a,alpha,max){
  Sm = 0
  m = 0
  continue = TRUE
  while(continue){
    m = m + 1
    k1 = qnorm(1-alpha/2,mean = 0,sd = 1)
    k2 = qnorm(alpha/2,mean = 0,sd = 1)
    
    ym1 = k1*sqrt(m)
    ym2 = k2*sqrt(m)
    
    Xm = rnorm(n = 1,mean = mu_a,sd = sigma_a)
    Sm = Sm + Xm
    
    if (Sm >= ym1 || Sm <= ym2 || m >= max){continue = FALSE}
  }
  return(c(M=m,Ym2=ym2,Zm=Sm,Ym1=ym1))
}

check.norm = function(mu_a,sigma_a,alpha,max){
  result = seq.norm(mu_a,sigma_a,alpha,max)
  reject = result[2] >= result[3] || result[3] >= result[4]
  return(reject)
}

B = 1000
tab.tp = array(data = NA, dim = c(length(mu_a),
                                  length(n),
                                  length(sigma_a),
                                  length(alpha)),dimnames = list(c(mu_a),
                                                                 c(n),
                                                                 c(sigma_a),
                                                                 c(alpha)))

for(l in 1:length(alpha)){
  for(k in 1:length(sigma_a)){
    for(j in 1:length(n)){
      for(i in 1:length(mu_a)){
        
        tab.tp[i,j,k,l] = sum(replicate(B,check.norm(alpha = alpha[l],
                                                     sigma_a = sigma_a[k],
                                                     max = n[j],
                                                     mu_a = mu_a[i])))/B
        
      }
    }
  }
}

tab.tp = tab.tp %>% melt(varnames = c("mu_a","n","sigma_a","alpha"))

labeller = label_bquote(rows = alpha ==. (alpha),
                        cols = sigma[a] ==. (sigma_a))

tab.tp %>% ggplot() +
  geom_line(mapping = aes(x = mu_a,y = value, col = as.factor(n))) +
  facet_grid(rows = vars(alpha), cols = vars(sigma_a),labeller = labeller) +
  theme_bw() + 
  labs(title = "Power Curves",
       subtitle = expression(paste("Power as a function of ",mu[a], ", ",
                                   sigma[a], ", Sample Size (n) and Significance Lvl (",alpha,") when ",
                                   mu[0]," = ", 0, " and ",sigma[0]," = ",1)),
                                   x = expression(paste(mu[a])),
                                   y = expression(paste("P[X" >= "x|",mu, "=", mu[a],"]")),
                                   color = "n")


# True Positive Rate - t-distribution ----
N = c(10,25,50,100,150) # vector of interim analyses
n = c(1,3,6) # vector of additional data points gathered 
alpha = c(.10,.05,.01) # vector of significance levels
mu_a = c(-2,-1,0,1,2) # vector of alternative means
sigma_a = c(.5,1,2,3,10) # vector of alternative sigmas

# N = 10 # of interim analyses
n = 5 # of data gathered for each interim analysis
# alpha = 0.05
# mu_a = 0
# sigma_a = 1

seq.t = function(mu_a,sigma_a,alpha,n,max){
  Sm = 0
  m = 0 # tracks the # of interim analyses
  k = 0 # tracks the total # of data points
  x_i = numeric()
  continue = TRUE
  while(continue){
    m = m + 1
    k = k + n
    #k1 = qnorm(1-alpha/2,mean = 0,sd = 1)
    #k2 = qnorm(alpha/2,mean = 0,sd = 1)
    
    #ym1 = k1*sqrt(m)
    #ym2 = k2*sqrt(m)
    
    Xm = rnorm(n = n,mean = mu_a,sd = sigma_a)
    x_i = append(x_i,Xm)
    Sm = Sm + sum(Xm)
    
    k1 = qt(p = 1-alpha/2,df = k-1)
    k2 = qt(p = alpha/2,df = k-1)

    ym1 = k1*sqrt(k)*sd(x_i)
    ym2 = k2*sqrt(k)*sd(x_i)
    
    if (Sm >= ym1 || Sm <= ym2 || m >= max){continue = FALSE}
  }
  return(c(M=m,K = k,Ym2=ym2,Zm=Sm,Ym1=ym1))
}

check.t = function(mu_a,sigma_a,alpha,n,max){
  result = seq.t(mu_a = mu_a,sigma_a = sigma_a,alpha,n = n,max = N)
  reject = result[3] >= result[4] || result[4] >= result[5]
  return(reject)
}

B = 1000
#sum(replicate(B,check.t(mu_a = mu_a,sigma_a = sigma_a,alpha = alpha,n = n,max = N)))/B

#tab.tp = matrix(data = NA,nrow = length(n),ncol = 1)

# for(i in 1:length(n)){
#   tab.tp[i] = sum(replicate(B,check.t(mu_a = mu_a,sigma_a = sigma_a,alpha = alpha,max = n)))/B
#   #tab.tp[i] = sum(replicate(B,check.norm(mu_a = mu_a,sigma_a = sigma_a,alpha = alpha,max = n)))/B
# }


tab.tp = array(data = NA, dim = c(length(mu_a),
                                  length(N),
                                  length(sigma_a),
                                  length(alpha)),dimnames = list(c(mu_a),
                                                                 c(N),
                                                                 c(sigma_a),
                                                                 c(alpha)))

for(l in 1:length(alpha)){
  for(k in 1:length(sigma_a)){
    for(j in 1:length(N)){
      for(i in 1:length(mu_a)){
        
        tab.tp[i,j,k,l] = sum(replicate(B,check.t(alpha = alpha[l],
                                                     sigma_a = sigma_a[k],
                                                     n = n,
                                                     max = N[j],
                                                     mu_a = mu_a[i])))/B
        
      }
    }
  }
}

tab.tp = tab.tp %>% melt(varnames = c("mu_a","N","sigma_a","alpha"))

labeller = label_bquote(rows = alpha ==. (alpha),
                        cols = sigma[a] ==. (sigma_a))

tab.tp %>% ggplot() +
  geom_line(mapping = aes(x = mu_a,y = value, col = as.factor(N))) +
  facet_grid(rows = vars(alpha), cols = vars(sigma_a),labeller = labeller) +
  theme_bw() + 
  labs(title = "Power Curves",
       subtitle = expression(paste("Power as a function of ",mu[a], ", ",
                                   sigma[a], ", Sample Size (n) and Significance Lvl (",alpha,") when ",
                                   mu[0]," = ", 0, " and ",sigma[0]," = ",1,
                                   " with using the t-distribution")),
       x = expression(paste(mu[a])),
       y = expression(paste("P[X" >= "x|",mu, "=", mu[a],"]")),
       color = "# of Interim Analyses")







# True Positive Rate - Regression ----

B0 = 10
B1 = 2
sigma = 1
nsim = 1000
pval = numeric(nsim)
N.vec = seq(25,100,by = 1)
power.N = numeric(length(N.vec))
for(j in 1:length(N.vec)){
  N = N.vec[j]
  x = seq(1,20,length = N.vec[j])
  for(i in 1:nsim){
    y_det = B0 + B1*x
    y = rnorm(N, mean = y_det, sd = sigma)
    m = lm(y ~ x)
    pval[i] = coef(summary(m))["x","Pr(>|t|)"]
  }
  power.N[j] = sum(pval < 0.05)/nsim
}
power.N
plot(N.vec,power.N,type = "l")

# aging surveillance study design simulation
# false positive rate
# true positive rate
# N interim analyses and n data points

# factors - alpha, sigma, B0, B1
# how to determine frequency of test? based on age
# need historical test sizes and ages

pval = numeric(nsim)
N.vec = c(10,25,50) # vector of interim analyses
power.N = numeric(length(N.vec)) # vector of power achieved at each # of interim analyses

power.sim = function(x.new,B0,B1){
    x.old = x_i
    y.old = y_i
    
    y_det = B0 + B1*x.new
    y.new = rnorm(n, mean = y_det, sd = sigma)
    
    x = append(x.old,x.new)
    y = append(y.old,y.new)
    
    m = lm(y ~ x)
    return(c(x = list(x), y = list(y), p.value = coef(summary(m))["x","Pr(>|t|)"]))
}

pop.N = 400
nsim = 1000
B0 = 10 # true intercept
B1 = 1 # true slope
sigma = 5 # true sd

t = 2 # test frequency in years
n = 4 # test quantity
power.N.vec = numeric(length(N.vec))  
for (j in 1:length(N.vec)){
  
  pop.age = seq(0,5,length.out = 1000) %>% round(digits = 2) %>% sample(pop.N)
  pop.age_indices = seq(1,pop.N,by = 1)
  
  N = N.vec[j]
  power.N = numeric(N)
  time = 0 # tracks the time
  m = 0 # tracks the running # of interim analyses
  k = 0 # tracks the running sample size
  x_i = numeric() # stores the data - age
  y_i = numeric() # stores the data - value
  continue = TRUE
  while(continue){
    m = m + 1
    k = k + n
    time = time + t
    pop.age = pop.age + t
    
    sampled.indices = sample(pop.age_indices,size = n,replace = FALSE) # sample the indices
    sample.x = pop.age[sampled.indices]
    
    power.N[m] = sum(replicate(n = nsim,expr = power.sim(x.new = sample.x,B0 = B0,B1 = B1)$p.value) < .01)/nsim
    
    x_i = append(x_i,sample.x)
    sample.y = rnorm(n, mean = B0 + B1*sample.x, sd = sigma)
    y_i = append(y_i,sample.y)
    pop.age_indices = pop.age_indices[-sampled.indices]
    
    if(m >= N){continue = FALSE}
  }
  power.N.vec[j] = power.N %>% tail(n = 1)
}  
  




  max = 20
  power.N = numeric(max)
  time = 0 # tracks the time
  m = 0 # tracks the running # of interim analyses
  k = 0 # tracks the running sample size
  x_i = numeric() # stores the data - age
  y_i = numeric() # stores the data - value
  continue = TRUE
  while(continue){
    m = m + 1
    k = k + n
    time = time + t
    pop.age = pop.age + t
    
    sampled.indices = sample(pop.age_indices,size = n,replace = FALSE) # sample the indices
    sample.x = pop.age[sampled.indices]
    
    power.N[m] = sum(replicate(n = nsim,expr = power.sim(x.new = sample.x,B0 = B0,B1 = B1)$p.value) < .01)/nsim
    
    x_i = append(x_i,sample.x)
    sample.y = rnorm(n, mean = B0 + B1*sample.x, sd = sigma)
    y_i = append(y_i,sample.y)
    pop.age_indices = pop.age_indices[-sampled.indices]
    
    if(m >= max){continue = FALSE}
  }
  
  
  
  

  
  
  power.N[j] = sum(pval < 0.05)/nsim
}

power.N
plot(N.vec,power.N,type = "l")
























