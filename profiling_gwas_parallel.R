
###################  function for GWAS  ##################

profiling_gwas <- function(P, Genotypes, K, ID, name="gwas"){
  start_time <- Sys.time()
  n = nrow(P)
  m = ncol(P)
  ######## eigen decomposition
  r <- eigen(K)
  U <- r$vectors
  r$values[r$values<1e-5] =1e-5
  D = r$values
  H = sqrt(solve(diag(D)))%*%t(U)
  
  ######### transformation
  
  Y <- H%*%P
  X<- H%*%rep(1,n)
  W = D
  Z <- H%*%t(Genotypes)
  
  data_wide = data.frame(
    ID = rep(1,n),
    X=X,
    W=W,
    Y=Y,
    Z=Z
  )
  # write.csv(data_wide, "data_wide_hybrid.csv",row.names = F)
  
  ID <- as.matrix(data_wide[, 1]) 
  X <- as.matrix(data_wide[,2])
  W = as.matrix(data_wide[,3])
  Y <- as.matrix(data_wide[, 4:(3+m)])
  Z0 <- as.matrix(data_wide[, (m+4):ncol(data_wide)])
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)
  print("Finish Eigen Decomposition Transformation")
  
  ######## EM algorithm to estimate aa and ee
  start_time <- Sys.time()
  aa <- diag(m)
  ee <- diag(m)
  maxiter <- 1000
  minerr <- 1e-8
  iter <- 0
  err <- 1e8  ## should update err
  www <- matrix(numeric(0), ncol=2+m)  # To store iteration results if needed
  
  while(iter < maxiter && err > minerr) {
    xx <- matrix(0, m, m)
    xy <- matrix(0, m, 1)
    yy <- matrix(0, m, m)
    bb <- matrix(0, m, m)
    
    for (j in 1:n) {
      Xj = diag(X[j],m)
      tt <- solve(aa * W[j] + ee+ diag(m) * 1e-5) #
      xx <- xx + W[j] * Xj %*% tt %*% Xj
      xy <- xy + W[j] * Xj %*% (tt %*% Y[j,])
    }
    
    mu <- solve(xx) %*% xy
    # another way to get mu (why different)
    # for (j in 1:n){
    #   xj = diag(X[j],m)
    #   tt = solve(aa+ (1/W[j])*ee)
    #   xx <- xx + Xj %*% tt %*% Xj
    #   xy <- xy + Xj %*% (tt %*% Y[j,])
    # }
    # mu <- solve(xx) %*% xy
    
    for (j in 1:n) {
      Xj = diag(X[j],m)
      tt <- solve(aa * W[j] + ee+ diag(m) * 1e-5 ) #+ diag(m) * 1e-5
      ex <- W[j] * aa %*% tt %*% (Y[j,] - Xj %*% mu)
      yy <- yy + W[j] * (Y[j,] - Xj %*% mu - ex)%*% Y[j,]
      va <- aa - W[j] * aa %*% tt %*% aa
      bb <- bb + (ex %*% t(ex) + va)
    }
    
    ee1 <- yy / (n - 1)
    aa1 <- bb / n
    
    ### info for each iteration
    www <- rbind(www, c(iter, err, t(mu)))#,ll
    err <- (sum((aa1 - aa)^2) + sum((ee1 - ee)^2)) / (2 * m^2)
    ee <- ee1
    aa <- aa1
    iter <- iter + 1
  }
  write.csv(aa,paste("G_EM_", name, ".csv",sep = ""),row.names = F)
  write.csv(ee,paste("R_EM_", name, ".csv",sep = ""),row.names = F)
  write.csv(mu,paste("beta_EM_", name, ".csv",sep = ""), row.names = F)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)
  print("Finish variance component estimation")
  
  
  ###### scan the genome 
  start_time <- Sys.time()
  lambda = (aa)%*%solve(ee)
  p <- ncol(Z0)
  n <- nrow(X)
  m <- ncol(Y)
  # result_scan <- matrix(numeric(0),ncol = m+4)  # To store iteration results if needed
  
  # parallel computing
  library(doParallel)
  library(foreach)
  
  # Detect number of available cores
  num_cores <- detectCores() - 14 # Leave one core free
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  p <- ncol(Z0)
  
  result_scan <- foreach(k = 1:p, .combine = rbind, .packages = "Matrix") %dopar% {
    
    Z <- Z0[,k]
    Axx <- matrix(0, m, m)
    Axz <- matrix(0, m, m)
    Axy <- matrix(0, m, 1)
    Azy <- matrix(0, m, 1)
    Azz <- matrix(0, m, m)
    Ayy <- matrix(0, m, m)
    
    for (j in 1:n) {
      Xj <- diag(X[j,], m)
      Zj <- diag(Z[j], m)
      tt <- solve(lambda * W[j] + diag(m) * (1 + 1e-5))
      
      Axx <- Axx + W[j] * Xj %*% tt %*% Xj
      Axy <- Axy + W[j] * Xj %*% (tt %*% (Y[j,]))
      Axz <- Axz + W[j] * Xj %*% tt %*% Zj
      Azy <- Azy + W[j] * Zj %*% (tt %*% (Y[j,]))
      Azz <- Azz + W[j] * Zj %*% tt %*% Zj
      Ayy <- Ayy + W[j] * (tt %*% (Y[j,])) %*% t(Y[j,])
    }
    
    Azx <- Axz
    Ayx <- Axy
    Czz <- Azz - Azx %*% solve(Axx) %*% Axz
    Czy <- Azy - Azx %*% solve(Axx) %*% Axy
    Cyy <- Ayy - (solve(Axx) %*% Axy) %*% t(Ayx)
    Cyz <- Czy
    
    gg <- solve(Czz) %*% Czy
    Crr <- Cyy - solve(Czz) %*% Czy %*% t(Cyz)
    ee <- Crr / (n - 2)
    vv <- solve(Czz) %*% ee
    vvi <- solve(vv)
    wald <- t(gg) %*% vvi %*% gg
    pvalue <- exp(pchisq(wald,df=m, lower.tail = FALSE, log.p = TRUE))
    log10p <- -log10(pvalue)
    
    # Return the results as a vector
    c(k, t(gg), wald, pvalue, log10p)
  }
  
  # Stop the parallel cluster
  stopCluster(cl)
  colnames(result_scan) <- c("snp", colnames(Y), "wald", "p", "log10p")
  
  result_scan = as.data.frame(result_scan)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)
  write.csv(result_scan,paste("results_profileGWAS_", name, ".csv",sep = ""),row.names = FALSE)
  return(result_scan)
}
