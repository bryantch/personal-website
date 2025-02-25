# power analysis for aging study design

library(ggplot2)
library(magrittr)
library(reshape2)

B0 = 20 # true intercept
B1 = c(.1,.5,1,2,5) # vector of true slopes
sigma = c(.1,.5,1,2,5) # vector of true sds
max = 10
t.vec = c(1,2,3) # vector test frequencies in years
n.vec = c(2,4,6,8) # vector of test quantities

slope.detect = function(x.new, B0, B1, sigma){
  x.old = x_i
  y.old = y_i
  
  y_det = B0 + B1*x.new
  y.new = rnorm(n = n, mean = y_det, sd = sigma)
  
  x = append(x.old, x.new)
  y = append(y.old, y.new)
  model = lm(y ~ x)
  return(c(x = list(x), y = list(y),p.value = coef(summary(model))["x","Pr(>|t|)"]))
}

power.sim = function(x.new,B0,B1,sigma,nsim){
  result = sum(replicate(n = nsim,expr = slope.detect(x.new = sample.x,
                                                      B0 = B0,
                                                      B1 = B1,
                                                      sigma = sigma)$p.value) < .05)/nsim
  return(result)
}

table.power = array(data = NA, dim = c(length(n.vec),
                                       length(t.vec),
                                       max,
                                       length(B1),
                                       length(sigma)),dimnames = list(c(n.vec),
                                                                      c(t.vec),
                                                                      c(seq(1,max,by = 1)),
                                                                      c(B1),
                                                                      c(sigma)))
nsim = 5000
for (l in 1:length(n.vec)) {
  for (c in 1:length(t.vec)) {
    for(i in 1:length(B1)){
      for(j in 1:length(sigma)){
        
        t = t.vec[c]
        n = n.vec[l]
        
        continue = TRUE
        m = 0
        k = 0
        time = 0
        x_i = numeric()
        y_i = numeric()
        
        pop.N = 400
        pop.age = seq(0,5,length.out = 1000) %>% sample(pop.N) %>% round(digits = 2)
        pop.age_indices = seq(1,pop.N,by = 1)
        
        while (continue) {
          m = m + 1
          k = k + n
          time = time + t
          pop.age = pop.age + t
          
          sampled.indices = sample(pop.age_indices, size = n, replace = FALSE) # sampling the indices
          sample.x = pop.age[sampled.indices]
          
          table.power[l,c,m,i,j] = power.sim(x.new = sample.x,B0 = B0,
                                             B1 = B1[i],
                                             sigma = sigma[j],
                                             nsim = nsim)
          
          x_i = append(x_i, sample.x)
          sample.y = rnorm(n = n,mean = B0 + B1[i]*sample.x,sd = sigma[j])
          y_i = append(y_i,sample.y)
          pop.age_indices = pop.age_indices[-sampled.indices]
          
          if(m >= max){continue = FALSE}
        }
        
      }
    }
  }
}

table.power 


labeller = label_bquote(cols = sigma ==. (sigma),
                        rows = beta[1] ==. (B1))

table.power %>% melt(varnames = c("n",
                                  "t",
                                  "N",
                                  "B1",
                                  "sigma")) %>% ggplot() +
  geom_line(mapping = aes(x = N,y = value, col = as.factor(n),linetype = as.factor(t))) +
  facet_grid(rows = vars(B1), cols = vars(sigma),labeller = labeller) +
  theme_bw() +
  labs(title = "Power Curves",
       subtitle = expression(paste("Power as a function of ",beta[1], ", ",
                                   sigma, ", Test Size (n), Test Frequency (t), # of Interim Analyses (N) and Significance Lvl (",alpha,")")),
       x = expression(paste(N)),
       y = expression(paste("P[T" >= "t|",beta[1], " > " , 0,"]")),
       color = "Test Quantity",linetype = "Test Frequency")






