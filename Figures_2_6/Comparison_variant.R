library("vioplot")
source("Functions_priori.R")

K <- 50
alpha <- c(1,rep(1/(K-1), K-1))
n <- 1000


#preprocessing
trials <- 300
iterations <- ceiling(50)
#store first component
P_nrstore1 <- matrix(0, trials, iterations)
P_nrstore_variant1 <- matrix(0, trials, iterations)
P_nrstore2 <- matrix(0, trials, iterations)
P_nrstore_variant2 <- matrix(0, trials, iterations)

for(t in 1:trials){
  
  #initial partition
  c <- sample(1:K, n, replace = T)
  #c <- rep(2, n)
  
  #non reversible
  res_nr <- nonrev_sampler_priori(alpha, n, K, n*iterations, c, verbose = F)
  S_nr <- res_nr[[1]]
  P_nrstore1[t,] <- S_nr[seq(1, n*iterations, floor(n^(1))),1]/n
  P_nrstore2[t,] <- S_nr[seq(1, n*iterations, floor(n^(1))),2]/n
  
  #variant
  res_nr_variant <- nonrev_sampler_priori_variant(alpha, n, K, n*iterations, c)
  S_nr_variant <- res_nr_variant[[1]]
  P_nrstore_variant1[t,] <- S_nr_variant[seq(1, n*iterations, floor(n^(1))),1]/n
  P_nrstore_variant2[t,] <- S_nr_variant[seq(1, n*iterations, floor(n^(1))),2]/n
  
  if(t %% 50 == 0){
    print(paste("Trial", t, "done!"))
    #print(round(to_print1, 2))
    #print(round(to_print2, 2))
    #cat("\n")
  }
  
}

#non reversible
plot(P_nrstore1[,iterations], P_nrstore2[,iterations], pch = 19, main = "Non reversible", ylab = "", xlab = "", xlim = c(0,1), ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
hist(P_nrstore1[,iterations], main = "Non reversible", freq = F, ylab = "", xlab = "", col = "red", xlim = c(0,1))
curve(dbeta(x,alpha[1], sum(alpha[-1])), col = "blue", lwd = 2, add = T)

#non reversible
plot(P_nrstore_variant1[,iterations], P_nrstore_variant2[,iterations], pch = 19, main = "Non reversible variant", ylab = "", xlab = "", xlim = c(0,1), ylim = c(0,1), cex.axis = 1.5, cex.main = 1.5)
hist(P_nrstore_variant1[,iterations], main = "Non reversible", freq = F, ylab = "", xlab = "", col = "red", xlim = c(0,1))
curve(dbeta(x,alpha[1], sum(alpha[-1])), col = "blue", lwd = 2, add = T)

#vioplot
ind <- c(1, seq(10, iterations, 5))
to_plot_nr <- P_nrstore1[,ind]
colnames(to_plot_nr) <- ind
to_plot_nr <- as.data.frame(to_plot_nr)
to_plot_nr_variant <- P_nrstore_variant1[,ind]
colnames(to_plot_nr_variant) <- ind
to_plot_nr_variant <- as.data.frame(to_plot_nr_variant)

vioplot(to_plot_nr, side = "left", col = "black", border = "NA",  drawRect = F, main = "Marginal distribution of the chains", cex.axis = 1.5, cex.main = 1.5)
vioplot(to_plot_nr_variant, side = "right", col = "gray", border = "NA",  drawRect = F, add = T)
  
