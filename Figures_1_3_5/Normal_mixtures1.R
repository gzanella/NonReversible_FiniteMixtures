source("https://raw.githubusercontent.com/gzanella/NonReversible_FiniteMixtures/refs/heads/main/Figures_1_3_5/Functions_normal.R")
library("vioplot")


#single example
K <- 2
alpha <- rep(0.5, K)
n <- 2000

#prior
mu0 <- 0
sigma20 <- 1
sigma2 <- 1

#generate data
w_true <- c(0.9, 0.1)
mu_true <- c(0.9, -0.9)
set.seed(1920)
y <- exact_sample_mixture(w_true, n, K, mu_true, rep(sigma2, K))
plot(density(y))
curve(dnorm(x, 0.8, sqrt(sigma2)), add = T, col = 'red')

#iterations
iterations <- 150 #ceiling(15*log(n))

trials <- 100
Store_r <- matrix(0, trials, iterations)
Store_nr <- matrix(0, trials, iterations)
for(i in 1:trials){
  c <- sample(1:K, n, replace = T)
  
  res_r <- rev_sampler(alpha, n, K, n*iterations, c, y, sigma2, mu0, sigma20)
  S_r <- res_r[[1]]
  Store_r[i, ] <- apply(S_r[seq(1, n*iterations, floor(n^(1))),], 1, max)
  
  res_nr <- nonrev_sampler(alpha, n, K, n*iterations, c, y, sigma2, mu0, sigma20, verbose = F)
  S_nr <- res_nr[[1]]
  Store_nr[i, ] <- apply(S_nr[seq(1, n*iterations, floor(n^(1))),], 1, max)
  
  if(i %% 5 == 0){
    print(paste("Trial", i, "done!"))
    #print(round(to_print1, 2))
    #print(round(to_print2, 2))
    #cat("\n")
  }
}
#save
#saveRDS(Store_r, file = "Example_normal_r.rds")
#saveRDS(Store_nr, file = "Example_normal_nr.rds")

#plot reversible
plot(Store_r[1, ]/n, type = "l", col = "black", lwd = 3, ylim = c(0.5,1), main = "Marginal", ylab = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
for(i in 2:trials){
  points(Store_r[i, ]/n, type = "l", col = gray(0.6), lwd = 3, ylim = c(0.5,1),  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
}
points(Store_r[1, ]/n, type = "l", col = "black", lwd = 3, ylim = c(0.5,1),  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
#plot non reversible
plot(Store_nr[1, ]/n, type = "l", col = "black", lwd = 3, ylim = c(0.5,1), main = "Non reversible", ylab = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
for(i in 2:trials){
  points(Store_nr[i, ]/n, type = "l", col = gray(0.5), lwd = 3, ylim = c(0.5,1),  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
}
points(Store_nr[1, ]/n, type = "l", col = "black", lwd = 3, ylim = c(0.5,1),  cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

#vioplot
ind <- c(1, seq(10, iterations, 10))
to_plot_r <- Store_r[, ind]/n
colnames(to_plot_r) <- ind
to_plot_r <- as.data.frame(to_plot_r)
to_plot_nr <- Store_nr[, ind]/n
colnames(to_plot_nr) <- ind
to_plot_nr <- as.data.frame(to_plot_nr)

vioplot(to_plot_nr, side = "left", col = "black", border = "NA",  drawRect = F, main = "Marginal distribution of the chains", cex.axis = 1.5, cex.main = 1.5)
vioplot(to_plot_r, side = "right", col = "gray", border = "NA",  drawRect = F, add = T)
