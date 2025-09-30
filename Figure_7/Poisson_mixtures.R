source("Functions_poisson.R")
library("vioplot")
#Large simulation generating from the Bayesian model

K <- 3
alpha <- rep(1, K)
n <- 1000

#prior
beta1 <- 1
beta2 <- 1


#preprocessing
trials <- 300
iterations <- ceiling(100)
#store first component
P_rstore1 <- matrix(0, trials, iterations)
P_nrstore1 <- matrix(0, trials, iterations)
P_rstore2 <- matrix(0, trials, iterations)
P_nrstore2 <- matrix(0, trials, iterations)
for(t in 1:trials){
  #generate data
  Lambda_true <- rgamma(K, beta1, beta2)
  p_true <- rgamma(K, alpha, 1)
  y <- exact_sample_mixture(p_true, n, K, Lambda_true)
  
  #initial partition
  #c <- 1 + rbinom(n, 1, 1/2)
  c <- sample(1:K, n, T)
  #c <- kmeans(y, centers = K)$cluster
  #c <- rep(1, n)
  
  #reversible
  res_r <- rev_sampler(alpha, n, K, n*iterations, c, y, beta1, beta2)
  S_r <- res_r[[1]]
  P_rstore1[t,] <- S_r[seq(1, n*iterations, floor(n^(1))),1]/n
  P_rstore2[t,] <- S_r[seq(1, n*iterations, floor(n^(1))),2]/n
  
  #non reversible
  res_nr <- nonrev_sampler(alpha, n, K, n*iterations, c, y, beta1, beta2, verbose = F)
  S_nr <- res_nr[[1]]
  P_nrstore1[t,] <- S_nr[seq(1, n*iterations, floor(n^(1))),1]/n
  P_nrstore2[t,] <- S_nr[seq(1, n*iterations, floor(n^(1))),2]/n
  
  if(t %% 50 == 0){
    print(paste("Trial", t, "done!"))
    #print(round(to_print1, 2))
    #print(round(to_print2, 2))
    #cat("\n")
  }
  
}
#save
saveRDS(P_rstore1, file = "Poisson_mixtures_alpha1_r1.rds")
saveRDS(P_rstore2, file = "Poisson_mixtures_alpha1_r2.rds")
saveRDS(P_nrstore1, file = "Poisson_mixtures_alpha1_nr1.rds")
saveRDS(P_nrstore2, file = "Poisson_mixtures_alpha1_nr2.rds")
#plots
#1)
#reversible
plot(P_rstore1[,iterations], P_rstore2[,iterations], pch = 19, main = "Marginal", ylab = "", xlab = "", xlim = c(0,1), ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
hist(P_rstore1[,iterations], main = "Reversible", freq = F, ylab = "", xlab = "", col = "red", xlim = c(0,1))
curve(dbeta(x,alpha[1], sum(alpha[-1])), col = "blue", lwd = 2, add = T)
#non reversible
plot(P_nrstore1[,iterations], P_nrstore2[,iterations], pch = 19, main = "Non reversible", ylab = "", xlab = "", xlim = c(0,1), ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
hist(P_nrstore1[,iterations], main = "Non reversible", freq = F, ylab = "", xlab = "", col = "red", xlim = c(0,1))
curve(dbeta(x,alpha[1], sum(alpha[-1])), col = "blue", lwd = 2, add = T)


#vioplot
ind <- c(1, seq(10, iterations, 10))
to_plot_r <- P_rstore1[,ind]
colnames(to_plot_r) <- ind
to_plot_r <- as.data.frame(to_plot_r)
to_plot_nr <- P_nrstore1[,ind]
colnames(to_plot_nr) <- ind
to_plot_nr <- as.data.frame(to_plot_nr)

vioplot(to_plot_nr, side = "left", col = "black", border = "NA",  drawRect = F, main = "Marginal distribution of the chains", cex.axis = 1.5, cex.main = 1.5)
vioplot(to_plot_r, side = "right", col = "gray", border = "NA",  drawRect = F, add = T)
