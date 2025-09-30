rev_sampler <- function(alpha, n, K, iterations, c, y, sigma2, mu0, sigma20){
  #K = number of clusters
  #alpha = K dimensional vector
  
  #preprocessing
  sums <- rep(0, K)
  tab <- table(c)
  sums[as.numeric(names(tab))] <- tab
  probs <- alpha/sum(alpha)
  
  sums_y <- rep(0, K)
  for(k in 1:K){
    sums_y[k] <- sum(y[c == k])
  }
  
  #to store
  Sums <- matrix(0, iterations, K)
  
  choices <- sample(1:n, iterations, replace = T)
  for(t in 1:iterations){
    
    i <- choices[t] #observation to update
    c_i <- c[i]
    y_i <- y[i]
    sums[c_i] <- sums[c_i] - 1
    sums_y[c_i] <- sums_y[c_i] - y_i
    
    #compute predictive means and variances
    means <- sigma20/(sigma2+sums*sigma20)*sums_y + sigma2/(sigma2+sums*sigma20)*mu0
    vars <- sigma2 + sigma20*sigma2/(sigma2+sums*sigma20)
    
    #sample new label
    log_prob <- log(alpha + sums) + dnorm(y_i, means, sqrt(vars), log = T)
    new_lab <- which.max(log_prob - log(-log(runif(K)))) #gumbel-max trick
    
    #update and store
    c[i] <- new_lab
    sums[new_lab] <- sums[new_lab] + 1
    sums_y[new_lab] <- sums_y[new_lab] + y_i
    
    Sums[t,] <- sums
    
  }
  
  return(list(Sums))
  
}
nonrev_sampler <- function(alpha, n, K, iterations, c, y, sigma2, mu0, sigma20, verbose = T){
  #K = number of clusters
  #alpha = K dimensional vector
  #sigma2 = variance of the likelihood
  #mu0 = mean of the baseline measure
  #sigma20 = variance of the baseline
  
  #preprocessing
  sums <- rep(0, K)
  tab <- table(c)
  sums[as.numeric(names(tab))] <- tab
  
  
  sums_y <- rep(0, K)
  for(k in 1:K){
    sums_y[k] <- sum(y[c == k])
  }
  #vector of velocities
  vel <- matrix(0, nrow = K, ncol = K)
  vel[lower.tri(vel)] <- sample(c(1, -1), K*(K-1)/2, replace = T)
  vel[upper.tri(vel)] <- - vel[lower.tri(vel)]
  
  #to store
  Sums <- matrix(0, iterations, K)
  #diagnostic
  counter_excursion <- 1
  excursions <- c()
  
  #random steps
  xi <- 1/2
  W <- runif(iterations) #for deciding whether to flip in the beginning
  to_flip_first <- log(W) < log(xi/n)
  U <- runif(iterations) #for deciding whether to accept in step 1
  V <- runif(iterations) #for deciding whether to flip in step 2
  to_flip <- log(V) < log(1-xi/n)
  V2 <- runif(iterations) #for choosing the observation in step 1
  
  #sample (k, k')
  U1 <- runif(iterations)
  U2 <- runif(iterations)
  i1 <- ceiling(n*U1) #to select the first component
  K_set <- 1:K
  i2 <- ceiling((K-1)*U2) #to select the second component
  
  for(t in 1:iterations){
    #sample clusters to consider
    k1 <- c[i1[t]]
    k2 <- (K_set[-k1])[i2[t]]
    
    #order the clusters
    velocity <- vel[k1, k2]
    velocity_old <- velocity
    #check if flip
    if(to_flip_first[t]){
      velocity <- - velocity
    }
    vel[k1, k2] <- velocity
    vel[k2, k1] <- - velocity
    #order the velocity
    if(velocity == 1){
      k_plus <- k1
      k_minus <- k2
    }
    else{
      k_plus <- k2
      k_minus <- k1
    }
    #k_plus <- ifelse(velocity == 1, k1, k2)
    #k_minus <- ifelse(velocity == 1, k2, k1)
    n_plus <- sums[k_plus]
    n_minus <- sums[k_minus]
    
    #1) propose a change
    if(n_minus > 0){
      #move one observation from k_minus to k_plus
      i <- ((1:n)[c == k_minus])[ceiling(n_minus*V2[t])]
      y_i <- y[i]
      
      #compute predictive mean and variance
      mean_plus <- sigma20/(sigma2+n_plus*sigma20)*sums_y[k_plus] + sigma2/(sigma2+n_plus*sigma20)*mu0
      var_plus <- sigma2 + sigma20*sigma2/(sigma2+n_plus*sigma20)
      mean_minus <- sigma20/(sigma2+(n_minus-1)*sigma20)*(sums_y[k_minus] - y_i) + sigma2/(sigma2+(n_minus-1)*sigma20)*mu0
      var_minus <- sigma2 + sigma20*sigma2/(sigma2+(n_minus-1)*sigma20)
      
      #prior
      log_alpha <- log(alpha[k_plus] + n_plus) + log(n_minus)-
        log(alpha[k_minus] + n_minus - 1) - log(n_plus + 1)
      
      #likelihood
      log_alpha <- log_alpha + 
        dnorm(y_i, mean_plus, sqrt(var_plus), log = T) - dnorm(y_i, mean_minus, sqrt(var_minus), log = T)
      
      
      if(log(U[t]) <= log_alpha){
        c[i] <- k_plus
        sums[k_minus] <- n_minus - 1
        sums[k_plus] <- n_plus + 1
        sums_y[k_minus] <- sums_y[k_minus] - y_i
        sums_y[k_plus] <- sums_y[k_plus] + y_i
        
        #flip
        velocity <- - velocity
      }
    }
    
    #2) make the flip
    if(to_flip[t]){
      velocity <- - velocity
    }
    vel[k1, k2] <- velocity
    vel[k2, k1] <- - velocity
    
    #store
    Sums[t,] <- sums
    
    #save diagnostic swaps
    if(velocity != velocity_old){
      excursions <- c(excursions, counter_excursion)
      counter_excursion <- 0
    }
    counter_excursion <- counter_excursion + 1
    
    #print(c(sums, velocity == velocity_old))
    
  }
  
  #excursions and diagnostic
  excursions <- c(excursions, counter_excursion)
  if(verbose){
    plot(excursions, type = 'p', pch = 19, col = "red", main = "Ordered excursions", ylab = "", xlab = "", log = "y")
    hist(log(excursions), col = "red", main = "Excursions", freq = F, ylab = "", xlab = "")
    print(summary(excursions))
    #print(paste("Proportion of excursions:", length(excursions)/iterations))
  }
  
  return(list(Sums))
  
}
exact_sample_mixture <- function(probs, n, K, mu, sigma2){
  #sample components
  probs <- probs/sum(probs)
  comp <- sample(1:K, n, replace = T, prob = probs)
  
  #sample observations
  y <- mu[comp] + sqrt(sigma2[comp])*rnorm(n)
  return(y)
}