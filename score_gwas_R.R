
scoretest_gwas <- function(P, Genotypes, K, ID, name){
  start_time <- Sys.time()
  n = nrow(P)
  m = ncol(P)
  p <- nrow(Genotypes)
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
  G = aa
  R = ee
  
  start_time <- Sys.time()
  
  pvalue_all = c()
  score_all = c()
  
  for (j in 1:p){
    # print(j)
    Zg = Z0[,j]
    U_all = matrix(0,nrow = m, ncol = 1)
    S_all = matrix(0, nrow=m, ncol=m*1)
    for (i in 1:n){
      # print(i)
      Xi = diag(X[i],m)
      Yi = Y[i,]
      Vi_inv = solve(G + (1/W[i])*R) # txt matrix
      Ai = Vi_inv%*%(Yi-Xi%*%mu)    # A tx1 matrix
      # obtain all Zg (all markers for that individual)
      Zgi_all <- diag(Zg[i],m)
      # Use lapply to multiply each matrix by A and collect results
      U_i <-  Zgi_all%*% Ai
      
      # calculate score function
      S_i =  Zgi_all%*%Vi_inv%*%Zgi_all
      
      # add the values to the matrix
      U_all = U_all+U_i
      S_all = S_all + S_i
    }
    score = t(U_all)%*% solve(S_all)%*%U_all
    pvalue =  exp(pchisq(score,df=m, lower.tail = FALSE, log.p = TRUE))
    score_all = c(score_all, score)
    pvalue_all = c(pvalue_all, pvalue)
  }
  result_score_gwas = data.frame(score = score_all, pvalue = pvalue_all)
  write.csv(result_score_gwas,
            paste("results_score_GWAS_", name, ".csv",sep = ""),row.names = F)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)
  print("Finish score test GWAS")
  return(result_score_gwas)
}

