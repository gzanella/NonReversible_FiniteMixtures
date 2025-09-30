source("Functions_priori.R")
#Large simulation generating from the prior

K <- 3
alpha <- rep(1, 3)
n <- 1000


#preprocessing
trials <- 300
iterations <- ceiling(100)
#store first component
P_rstore1 <- matrix(0, trials, iterations)
P_nrstore1 <- matrix(0, trials, iterations)
P_rstore2 <- matrix(0, trials, iterations)
P_nrstore2 <- matrix(0, trials, iterations)
for(t in 1:trials){
  
  #initial partition
  c <- sample(1:K, n, replace = T)
  #c <- rep(1, n)
  
  #reversible
  res_r <- rev_sampler_priori(alpha, n, K, n*iterations, c)
  S_r <- res_r[[1]]
  P_rstore1[t,] <- S_r[seq(1, n*iterations, floor(n^(1))),1]/n
  P_rstore2[t,] <- S_r[seq(1, n*iterations, floor(n^(1))),2]/n
  
  #non reversible
  res_nr <- nonrev_sampler_priori(alpha, n, K, n*iterations, c, verbose = F)
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
#saveRDS(P_rstore1, file = "Prior_mixtures_alpha1_r1.rds")
#saveRDS(P_rstore2, file = "Prior_mixtures_alpha1_r2.rds")
#saveRDS(P_nrstore1, file = "Prior_mixtures_alpha1_nr1.rds")
#saveRDS(P_nrstore2, file = "Prior_mixtures_alpha1_nr2.rds")
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

#2)
plot(1:iterations, P_nrstore1[1,], lwd = 3, type = "l", col = "black", main = "Traceplots", ylab = "", xlab = "", ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
points(1:iterations, P_rstore1[1,], lwd = 3, type = "l", col = gray(0.7), main = "Traceplots", ylab = "", xlab = "", ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
for(t in 2:100){
  points(1:iterations, P_nrstore1[t,], lwd = 3, type = "l", col = "black", main = "Traceplots", ylab = "", xlab = "", ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
}
for(t in 2:100){
  points(1:iterations, P_rstore2[t,], lwd = 3, type = "l", col = gray(0.7), main = "Traceplots", ylab = "", xlab = "", ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
}

#vioplot
library("vioplot")
ind <- c(1, seq(10, iterations, 10))
to_plot_r <- P_rstore1[,ind]
colnames(to_plot_r) <- ind
to_plot_r <- as.data.frame(to_plot_r)
to_plot_nr <- P_nrstore1[,ind]
colnames(to_plot_nr) <- ind
to_plot_nr <- as.data.frame(to_plot_nr)

vioplot(to_plot_nr, side = "left", col = "black", border = "NA",  drawRect = F, main = "Marginal distribution of the chains", cex.axis = 1.5, cex.main = 1.5)
vioplot(to_plot_r, side = "right", col = "gray", border = "NA",  drawRect = F, add = T)
