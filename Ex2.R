#needs at least 16GB of RAM, plus around 1Gb for each core above 4.
#takes about 4 hours total to run everything
#Should scale well up to around 50-100 cores, but only reproducable
#when using clustersize = 10 as here, otherwise set to half the number of
#logical cores

#make sure you have rstan version 2.26 or later (this is later than the current 
#release version) as the stan syntax for arrays has changed and the code will 
#not run otherwise. Find instructions here:
#https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

library(MASS)
library(expm)
library(rstan)
library(parallel)
library(rstudioapi)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
RNGkind("L'Ecuyer-CMRG")

#unrounder to avoid rounding errors
unrounder <- function(x){
  if (x<10^-16){
    return(10^-16)
  }
  if (x>1-10^-16){
    return(1-10^-16)
  }
  return(x)
} 

#Data generator
genner <- function(I,J,R,sigma_a,sigma_b,rho,bias){
  #generate alpha and beta 
  sh1 <- diag(x=c(sigma_a^2,sigma_b^2),nrow=2,ncol=2)
  sh2 <-  matrix(c(0,rho*sigma_a*sigma_b,rho*sigma_a*sigma_b,0),nrow=2,ncol=2)
  Sigma <- sh1 + sh2
  temp <- mvrnorm(n=J,mu=rep(0,2),Sigma=Sigma)
  alpha <- exp(temp[,1])
  beta <- exp(temp[,2])
  remove(temp)
  params = list("I"=I,"J"=J,"R"=R,"sigma_a"=sigma_a,"sigma_b"=sigma_b,"rho"=rho,
              "bias"=bias,"alpha"=alpha,"beta"=beta)
  #generate Y
  y <- vapply(1:J,FUN.VALUE = rep(1,I),FUN = function(j){
    rbeta(n=I,shape1=alpha[j],shape2=beta[j])
  })
   
  #add bias (if req)
  if(!bias==0){ 
    for (j in 1:J){
      for (i in 1:I){
        if(runif(1,0,1)<bias){
          y[i,j] <- runif(1,0,1)
        }
      }
    } 
  }
  #prevent rounding to 0 or 1
  for (j in 1:J){
    for (i in 1:I){
      y[i,j]<- unrounder(y[i,j])
    }
  }
  #generate holdout data
  h <- vapply(1:J,FUN.VALUE = rep(1,R),FUN = function(j){
    rbeta(n=R,shape1=alpha[j],shape2=beta[j])
  })
  if(!bias==0){ 
    for (j in 1:J){
      for (r in 1:R){
        if(runif(1,0,1)<bias){
          h[r,j] <- runif(1,0,1)
        }
      }
    } 
  }
  for (j in 1:J){
    for (r in 1:R){
      h[r,j]<- unrounder(h[r,j])
    }
  }
  
  return(list(y,h,params))
}

#use stan to generate posterior draws
stanner2 <-function(model,stan_input){
  sampling(model,data=stan_input,chains=chains,iter=2000,save_warmup = F, 
           control = list("max_treedepth"=15,"adapt_delta"=0.999),
           pars = c("alpha","beta","sigma_a","sigma_b","rho","bias"))
  
}

#calculate elpd using posterior draws and holdout data
elpder2 <- function(fit,y_holdout,incl_bias=F,cl=NA,calc_v=F,return_lpds=F){
  if(class(fit)[[1]] == "stanfit"){
    draws <- rstan::extract(fit, pars = c("alpha","beta","bias"))
  }else{
    draws <- fit
  }
  I = length(y_holdout[,1])
  J = length(y_holdout[1,])
  S = dim(draws[[1]])[1]
  input <- lapply(1:S,FUN=function(s){
    alpha <- draws[[1]][s,]
    beta <- draws[[2]][s,]
    bias <- draws[[3]][s]
    return(list(alpha,beta,bias))
  })
  if(anyNA(cl)){
    log_pd = lapply(X=input,FUN = function(x){
      alpha <- x[[1]]
      beta <- x[[2]]
      bias <- x[[3]]
      temp = vector(length=I*J)
      for (j in 1:J){
        temp[((j-1)*I+1):(j*I)] = dbeta(x=y_holdout[,j],shape1 = alpha[j],
                                        shape2 = beta[j],log = T)
      }
      if(incl_bias){
        temp = log(exp(temp)*(1-bias)+bias)
      }
      return(temp)
    })
  }else{
    log_pd = parLapply(cl,X=input,fun = function(x){
      alpha <- x[[1]]
      beta <- x[[2]]
      bias <- x[[3]]
      temp = vector(length=I*J)
      for (j in 1:J){
        temp[((j-1)*I+1):(j*I)] = dbeta(x=y_holdout[,j],shape1 = alpha[j],
                                        shape2 = beta[j],log = T)
      }
      if(incl_bias){
        temp = log(exp(temp)*(1-bias)+bias)
      }
      return(temp)
    })
  }
  log_pd = do.call(rbind,log_pd)
  v = var(matrixStats::rowLogSumExps(log_pd) - log(ncol(log_pd)))
  lpds <- matrixStats::colLogSumExps(log_pd) - log(nrow(log_pd))
  elpd <- mean(lpds)
  se <- sd(lpds)/sqrt(length(lpds))
  if(return_lpds&calc_v){
    return(list(c(elpd,se),lpds,v))
  }
  if(calc_v){
    return(list(c(elpd,se),v))
  }
  if(return_lpds){
    return(list(c(elpd,se),lpds))
  }
  return(c(elpd,se))
}

elpd_comparer2 <- function(fit1,fit2,bias1=F,bias2=F,y_holdout){
  lpds1 <- elpder2(fit=fit1,y_holdout=y_holdout,incl_bias=bias1,return_lpds=T)[[2]]
  lpds2 <- elpder2(fit=fit2,y_holdout=y_holdout,incl_bias=bias2,return_lpds=T)[[2]]
  lpds_comp <- lpds1 - lpds2
  elpd_comp <- mean(lpds_comp)
  se_comp <- sd(lpds_comp)/sqrt(length(lpds_comp))
  return(c(elpd_comp,se_comp))
}

#function needed for running posterior bootstrap
logliker <- function(params,y){
  J = (length(params)-2)/2
  alpha = params[1:J]
  beta = params[(J+1):(2*J)]
  if(any(alpha<0)|any(beta<0)){
    return(-Inf)
  }
  sig_a = params[2*J+1]
  sig_b = params[2*J+2]
  if(!((0<=sig_a&sig_a<=sigma_a_upper)&(0<=sig_b&sig_b<=sigma_b_upper))){
    return(-Inf)
  }
  loglik = sum(sapply(1:J, FUN=function(j){
    logf <- sum(dbeta(y[,j],shape1=alpha[j],shape2=beta[j],log=T))
    logg <- dlnorm(alpha,sdlog=sig_a,log=T)+dlnorm(alpha,sdlog=sig_a,log=T)
    return(logf+logg)
  }))
  return(loglik)
}


#function needed for running posterior bootstrap
sumlogfer <- function(theta,y_slice,w,w_0,lambda){
  alpha = theta[1]
  beta = theta[2]
  sigma_a = lambda[1]
  sigma_b = lambda[2]
  if(sigma_a<=0|sigma_b<=0|alpha<=0|beta<=0){
    return(-Inf)
  }
  temp1 <- sapply(y_slice,dbeta,shape1=alpha,shape2=beta,log=T)
  temp2 <- c(dlnorm(alpha,sdlog=sigma_a,log=T),dlnorm(beta,sdlog=sigma_b,log=T))
  return(w%*%temp1+w_0%*%temp2)
}

#function needed for running posterior bootstrap
sumlogger <- function(lambda,theta_slice,v,v_0){
  sigma_a = lambda[1]
  sigma_b = lambda[2]
  if(sigma_a<=0|sigma_b<=0){
    return(-Inf)
  }
  temp1 <- sapply(1:J,FUN=function(j){
    alpha = theta_slice[j,1]
    beta = theta_slice[j,2]
    if(alpha<=0|beta<=0){
      return(-Inf)
    }
    return(dlnorm(alpha,sdlog=sigma_a,log=T)+dlnorm(beta,sdlog=sigma_b,log=T))
  })
  temp2 <- c(dunif(sigma_a,min=0,max=sigma_a_upper,log=T), 
             dunif(sigma_b,min=0,max=sigma_b_upper,log=T))
  return(v%*%temp1+v_0%*%temp2)
}

#extract diagnostics from stanfit
diagnoser <- function(fit){
  divergences <- get_num_divergent(fit)
  max_treedepth <- get_num_max_treedepth(fit)
  bfmi <- get_bfmi(fit)
  max_rhat <- max(summary(fit)$summary[,10])
  min_neff <- min(summary(fit)$summary[,9])
  output <- list("divergences"=divergences,"max_treedepth"=max_treedepth,
                 "bfmi"=bfmi,"max_rhat"=max_rhat,"min_neff"=min_neff)
  return(output)
} 

#function to generate data with given params and run tests
totaliser2 <- function(config){
  seed = config[[1]]
  sigma_a = config[[2]]
  sigma_b = config[[3]]
  rho = config[[4]]
  bias = config[[5]]
  #generate data
  set.seed(seed[1])
  data = genner(I=I,J=J,R=R,sigma_a=sigma_a,sigma_b=sigma_b,rho=rho,bias=bias)
  y = data[[1]]
  y_holdout = data[[2]]
  model <- stan_model(file="ex2 combo.stan")
  std_input <- list("I"=I,"J"=J,"sigma_a_upper"=sigma_a_upper,
                    "sigma_b_upper"=sigma_b_upper,"rho_prior"=0,"bias_prior"=0,
                    "eta_a"=1,"eta_b"=1,"zeta"=1,"y"=y)
  #setup cluster
  cl <- makeCluster(cluster_size)
  clusterExport(cl,varlist = list("stanner2","elpder2","sumlogfer","sumlogger",
                                  "I","J","sigma_a_upper","sigma_b_upper",
                                  "B","grid_scale","coarse_scale","postboot_N",
                                  "chains","logliker"))
  clusterExport(cl,varlist = list("y","y_holdout","model","std_input"),
                envir=environment())
  clusterEvalQ(cl, {
    library(rstan)
    library(rstudioapi)
    options(mc.cores = chains)
  })
  
  
  #permutations of standard bayes
  clusterSetRNGStream(cl,iseed=seed[4])
  if(skip_Bayes){
    print("Skipping Bayes")
    results = lapply(1:4,FUN=function(x){
      return(list(NA,NA))
    })
  }else{
    print("Starting Bayes")
    index = rep(list(std_input),4)
    index[[2]][[5]] <- 1 #include rho in model
    index[[3]][[6]] <- 1 #include bias in model
    index[[4]][[5]] <- 1 #include both rho
    index[[4]][[6]] <- 1 #and bias in model
    
    results <- parLapply(cl,X=index,fun=function(input){
      fit <- stanner2(model,input)
      elpd <- elpder2(fit,y_holdout,incl_bias=input[[6]])
      return(list(elpd,fit))
    })
    remove(index)
    print("Finished Bayes")
  }
  
  #Bayesbag M = I
  clusterSetRNGStream(cl,iseed=seed[5])
  if(B==0){ #skip Bayesbag if B=0
    print("Skipping M=I Bayesbag")
    bayesbag_fit <- NA
    bayesbag_elpd <- NA
  }else{
    print("Starting M = I Bayesbag")
    bayesbag_sflist <- parLapply(cl,(1:B)*4,fun = function(cid){
      y_bagged <- matrix(nrow=I,ncol=J)
      for (j in 1:J){
          y_bagged[,j] = sample(x=y[,j],size=I,replace=T)
      }
      input <- std_input
      input[[10]] <- y_bagged
      return(stanner2(model,input))
    })
    
    bayesbag_fit <- sflist2stanfit(bayesbag_sflist)
    bayesbag_elpd <- elpder2(bayesbag_fit,y_holdout,cl=cl)
    remove(bayesbag_sflist)
    print("Finished M = I Bayesbag")
  }
  results[[5]] <- list(bayesbag_elpd,bayesbag_fit)
  remove("bayesbag_elpd","bayesbag_fit")
  gc()
  
  #Bayesbag opt M
  clusterSetRNGStream(cl,iseed=seed[6])
  if(B==0|skip_Bayes){ #skip opt M Bayesbag if B=0 or skipped bayes
    print("Skipping opt M Bayesbag")
    bayesbag_opt_M <- NA
    bayesbag_opt_fit <- NA
    bayesbag_opt_elpd <- NA
  }else{
    print("Starting opt M Bayesbag")
    v_bay=elpder2(fit=results[[1]][[2]],y_holdout=y_holdout,calc_v=T,cl=cl)[[2]]
    v_bag=elpder2(fit=results[[5]][[2]],y_holdout=y_holdout,calc_v=T,cl=cl)[[2]]
    bayesbag_opt_M <- round(I*v_bag/(v_bag-v_bay))
    print(bayesbag_opt_M)
    remove(v_bay,v_bag)
    if(!bayesbag_opt_M>=I){
      print("Optimal M below I, so skipping opt M Bayesbag")
      bayesbag_opt_fit <- NA
      bayesbag_opt_elpd <- NA
    }else{
      if(bayesbag_opt_M>10*I){
        bayesbag_opt_M = 10*I
      }
      bayesbag_opt_sflist <- parLapply(cl,(1:B)*4,fun = function(cid){
        y_bagged <- matrix(nrow=bayesbag_opt_M,ncol=J)
        for (j in 1:J){
          y_bagged[,j] = sample(x=y[,j],size=bayesbag_opt_M,replace=T)
        }
        input <- std_input
        input[[1]] = bayesbag_opt_M
        input[[10]] <- y_bagged
        return(stanner2(model,input))
      })
      
      bayesbag_opt_fit <- sflist2stanfit(bayesbag_opt_sflist)
      bayesbag_opt_elpd <- elpder2(bayesbag_opt_fit,y_holdout,cl=cl)
      remove(bayesbag_opt_sflist)
      print("Finished Opt M Bayesbag")
    }
  }
  results[[6]] <- list(bayesbag_opt_elpd,bayesbag_opt_fit,bayesbag_opt_M)
  remove("bayesbag_opt_elpd","bayesbag_opt_fit","bayesbag_opt_M")
  gc()
  
  #SMI
  clusterSetRNGStream(cl,iseed=seed[7])
  if(grid_scale==0){ #skip SMI if grid_scale=0
    print("Skipping SMI")
    SMI_fit <- NA
    SMI_elpd <- NA
    SMI_eta <- NA
    SMI_all_elpds <- NA
  }else{
    print("Starting SMI")
    grid <- list()
    for (p in 0:(grid_scale-1)){
      for (q in 0:(grid_scale-1)){
        grid[[p*grid_scale+q+1]] = c(p,q)
      }
    }
    holder <- parLapply(cl,grid,fun=function(coords){
      gc()
      p = coords[1]
      q = coords[2]
      eta_a=p/(grid_scale-1)
      eta_b=q/(grid_scale-1)
      input <- std_input
      input[[7]] <- eta_a
      input[[8]] <- eta_b
      fit <- stanner2(model,input)
      temp <- elpder2(fit,y_holdout,return_lpds = T)
      elpd <- temp[[1]]
      lpds <- temp[[2]]
      return(list(elpd,fit,lpds))
    })
    
    SMI_all_elpds <- array(dim=c(grid_scale,grid_scale,2))
    for (p in 0:(grid_scale-1)){
      for (q in 0:(grid_scale-1)){
        SMI_all_elpds[p+1,q+1,] = holder[[p*grid_scale+q+1]][[1]]
      }
    }
    
    SMI_opt_coords = which(SMI_all_elpds[,,1]==max(SMI_all_elpds[,,1]),arr.ind=T)
    trans_coords = (SMI_opt_coords[1]-1)*grid_scale+SMI_opt_coords[2]
    t1 <- holder[[trans_coords]][[3]]
    t2 <- holder[[grid_scale^2]][[3]]
    t3 = t1 - t2
    t4 = mean(t3) - 1.96*sd(t3)/sqrt(length(t3))
    if(t4<=0){
      SMI_opt_coords = c(grid_scale,grid_scale)
      trans_coords = (SMI_opt_coords[1]-1)*grid_scale+SMI_opt_coords[2]
    }
    SMI_elpd <- holder[[trans_coords]][[1]]
    SMI_fit <- holder[[trans_coords]][[2]]
    SMI_eta <- (SMI_opt_coords-1)/(grid_scale-1)
    remove("holder","SMI_opt_coords","trans_coords","grid","p","q",
           "t1","t2","t3","t4")
    print("Finished SMI")
  }
  results[[7]] <- list(SMI_elpd,SMI_fit,SMI_eta,SMI_all_elpds)
  remove("SMI_elpd","SMI_fit","SMI_eta","SMI_all_elpds")
  gc()
  
  #Coarsening
  clusterSetRNGStream(cl,iseed=seed[8])
  if(grid_scale==0){ #skip Coarsening if coarse_scale=0
    print("Skipping Coarsening")
    coarse_fit <- NA
    coarse_elpd <- NA
    coarse_zeta <- NA
    coarse_all_elpds <- NA
  }else{
    print("Starting Coarsening")
    holder <- parLapply(cl,(0:(coarse_scale-1)),fun=function(x){
      zeta <- x/(coarse_scale-1)
      input <- std_input
      input[[9]] = zeta
      fit <- stanner2(model,input)
      temp <- elpder2(fit,y_holdout,return_lpds = T)
      elpd = temp[[1]]
      lpds = temp[[2]]
      return(list(elpd,fit,lpds))
    })
    
    coarse_all_elpds <- matrix(ncol=2,nrow=coarse_scale)
    for (x in 1:coarse_scale){
      coarse_all_elpds[x,] = holder[[x]][[1]]
    }
    coarse_opt_coord=which(coarse_all_elpds[,1]==max(coarse_all_elpds[,1]),arr.ind=T)
    t1<-holder[[coarse_opt_coord]][[3]]
    t2<-holder[[coarse_scale]][[3]]
    t3 <- t1 - t2
    t4 <- mean(t3) - 1.96*sd(t3)/sqrt(length(t3))
    if(t4 <= 0){
      coarse_opt_coord <- coarse_scale
    }
    coarse_elpd <- holder[[coarse_opt_coord]][[1]]
    coarse_fit <- holder[[coarse_opt_coord]][[2]]
    coarse_zeta <- (coarse_opt_coord-1)/(coarse_scale-1)
    remove(holder,coarse_opt_coord,t1,t2,t3,t4)
    print("Finished Coarsening")
  }
  results[[8]] <- list(coarse_elpd,coarse_fit,coarse_zeta,coarse_all_elpds)
  remove("coarse_elpd","coarse_fit","coarse_zeta","coarse_all_elpds")
  gc()
  
  #Postboot
  clusterSetRNGStream(cl,iseed=seed[9])
  if(postboot_N==0){ #skip Posterior Bootstrap if postboot_N=0
    print("Skipping Posterior Bootstrap")
    postboot_fit <- NA
    postboot_elpd <- NA
    postboot_weights <- NA
    postboot_MLEs <- NA
  }else{
    print("Starting Posterior Bootstrap")
    #find mles
    MLEs <- optim(par=rep(1,2*J+2),fn=logliker,y=y,
                  control=list("fnscale"=-1,"maxit"=10^9))[[1]]
    theta_MLE <- matrix(nrow=J,ncol=2)
    theta_MLE[,1] = MLEs[1:J]
    theta_MLE[,2] = MLEs[(J+1):(2*J)]
    lambda_MLE = MLEs[(2*J+1):(2*J+2)]

    #use mles to estimate fisher information, use this to calc optimal weights
    w_0 <- t(vapply(X=1:J,FUN.VALUE=c(1,1),FUN=function(j){
      alpha = theta_MLE[j,1]
      beta = theta_MLE[j,2]
      I_f = Reduce('+',lapply(X=1:I,FUN=function(i){
        I_f_temp = matrix(ncol=2,nrow=2)
        t1 =log(y[i,j])-digamma(alpha)+digamma(alpha+beta)
        t2 = log(1-y[i,j])-digamma(beta)+digamma(alpha+beta)
        I_f_temp[1,1] = t1^2
        I_f_temp[1,2] = I_f_temp[2,1] = t1*t2
        I_f_temp[2,2] = t2^2
        return(I_f_temp)
      }))/I
      sqrt_I_f = sqrtm(I_f)
      
      J_f =  matrix(-trigamma(alpha+beta),ncol=2,nrow=2)
      J_f[1,1] = J_f[1,1] + trigamma(alpha)
      J_f[2,2] = J_f[2,2] + trigamma(beta)
      inv_J_f = solve(J_f)
      
      return(diag(sqrt_I_f %*% inv_J_f %*% sqrt_I_f))
    }))
    
    sig_a = lambda_MLE[1]
    sig_b = lambda_MLE[2]
    I_g = Reduce('+',lapply(X=1:J,FUN=function(j){
      alpha = theta_MLE[j,1]
      beta = theta_MLE[j,2]
      t1 = ((log(alpha)/sig_a)^2-1)/sig_a
      t2 = ((log(beta)/sig_b)^2-1)/sig_b
      I_g_temp = matrix(ncol=2,nrow=2)
      I_g_temp[1,1] = t1^2
      I_g_temp[1,2] = I_g_temp[2,1] = t1*t2
      I_g_temp[2,2] = t2^2
      return(I_g_temp)
    }))/J
    sqrt_I_g = sqrtm(I_g)
    
    J_g = Reduce('+',lapply(X=1:J,FUN=function(j){
      alpha = theta_MLE[j,1]
      beta = theta_MLE[j,2]
      J_g_temp = matrix(0,ncol=2,nrow=2)
      J_g_temp[1,1] = (3*(log(alpha)/sig_a)^2-1)/(sig_a^2)
      J_g_temp[2,2] = (3*(log(beta)/sig_b)^2-1)/(sig_b^2)
      return(J_g_temp)
    }))/J
    inv_J_g = diag(x=c(1/J_g[1,1],1/J_g[2,2]))
    v_0 <-diag(sqrt_I_g %*% inv_J_g %*% sqrt_I_g)
    remove(MLEs,lambda_MLE,J_g,inv_J_g,I_g,sqrt_I_g,sig_a,sig_b)
    
    #run double postboot alorithm using optimal weights
    
    #find data-liklihood-only mle for alpha and beta
    theta_hat <- parLapply(cl,1:J,fun=function(j){
      optim(par=c(1,1),fn=sumlogfer,y_slice=y[,j],w=rep(1,I),w_0=c(0,0),
            lambda=c(1,1),control=list("fnscale"=-1,"maxit"=10^9))[[1]]
    })
    theta_hat = do.call(rbind,theta_hat)
    
    clusterExport(cl,varlist=list("theta_hat","w_0","v_0"),envir = environment())
    
    
    draws <- parLapply(cl,1:postboot_N,fun=function(k){
      lambda_hat <- optim(par=c(1,1),fn=sumlogger,theta_slice=theta_hat,v=rexp(J,1),
                          v_0=v_0,control=list("fnscale"=-1,"maxit"=10^9))[[1]]
      theta <- t(vapply(1:J,FUN.VALUE = c(1,1),FUN = function(j){
        w = rexp(I,1)
        optim(par=c(1,1),fn=sumlogfer,y_slice=y[,j],w=w,w_0=w_0[j,],
              lambda=lambda_hat,control=list("fnscale"=-1,"maxit"=10^9))[[1]]
      }))
      lambda <- optim(par=c(1,1),fn=sumlogger,theta_slice=theta,v=rexp(J,1),
                      v_0=v_0,control=list("fnscale"=-1,"maxit"=10^9))[[1]]
      return(c(lambda,theta[,1],theta[,2]))
    })
    draws=do.call(rbind,draws)
    postboot_fit <- list("alpha" = draws[,3:12],"beta" = draws[,13:22],
                         "sigma_a" = draws[,1],"sigma_b" = draws[,2])
    postboot_elpd <- elpder2(postboot_fit,y_holdout)
    postboot_weights <- list("w_0"=w_0,"v_0"=v_0)
    postboot_MLEs <- list("theta_MLE"=theta_MLE,"theta_hat" = theta_hat)
    remove("theta_hat","theta_MLE","draws","w_0","v_0")
    print("Finished Posterior Bootstrap")
  }
  results[[9]] <-list(postboot_elpd,postboot_fit,postboot_weights,postboot_MLEs)
  remove("postboot_elpd","postboot_fit","postboot_weights")
  gc()
  
  stopCluster(cl)
  
  names(results) <- c("Bayes","corr Bayes","bias Bayes", "corr-bias Bayes",
                      "Bayesbag","SMI","Coarsening","Posterior Boostrap")
  return(list(results,data))
}

#set parameters (shared by all)
R=1000
I = 5
J = 10
sigma_a_upper = 5
sigma_b_upper = 5
skip_Bayes = F
B = 50
grid_scale = 11
coarse_scale=11
postboot_N = 2000
cluster_size = 10
chains = 2

#gen seeds
set.seed(7395723)
seeds = lapply(1:10,FUN=function(x){
  return(sample(1:10^9,10,replace=T))
})

#create specifications for each run (vary bias, corr etc.)
config <- lapply(1:4,FUN = function(x){
  list("seed"=seeds[[x]],"sigma_a"=3/2,"sigma_b"=3/2,"rho"=0,"bias"=0)
})
config[[2]][[4]] = 0.99
config[[3]][[5]] = 0.75
config[[4]][[4]] = 0.99
config[[4]][[5]] = 0.75

results <- lapply(config,totaliser2)
save(list="results",file="ex2_results.RData")

#prints diagnostic information
lapply(1:8,FUN=function(x){
 diagnoser(results[[1]][[1]][[x]][[2]])
})
lapply(1:8,FUN=function(x){
  diagnoser(results[[2]][[1]][[x]][[2]])
})
lapply(1:8,FUN=function(x){
  diagnoser(results[[3]][[1]][[x]][[2]])
})
lapply(1:8,FUN=function(x){
  diagnoser(results[[4]][[1]][[x]][[2]])
})


#section to create graphs
library(ggplot2)
library(tidyverse)
library(cowplot)
library(Cairo)

good_elpds <- lapply(results[[1]][[1]],FUN=function(x){
  return(x[[1]])
})
corr_elpds <- lapply(results[[2]][[1]],FUN=function(x){
  return(x[[1]])
})
bias_elpds <- lapply(results[[3]][[1]],FUN=function(x){
  return(x[[1]])
})
corr_bias_elpds <- lapply(results[[4]][[1]],FUN=function(x){
  return(x[[1]])
})

elpds_list <- list(good_elpds,corr_elpds,bias_elpds,corr_bias_elpds)

good_fits <- lapply(results[[1]][[1]],FUN=function(x){
  return(x[[2]])
})
corr_fits <- lapply(results[[2]][[1]],FUN=function(x){
  return(x[[2]])
})
bias_fits <- lapply(results[[3]][[1]],FUN=function(x){
  return(x[[2]])
})
corr_bias_fits <- lapply(results[[4]][[1]],FUN=function(x){
  return(x[[2]])
})

fits_list <- list(good_fits,corr_fits,bias_fits,corr_bias_fits)

good_y_holdout <- results[[1]][[2]][[2]]
corr_y_holdout <- results[[2]][[2]][[2]]
bias_y_holdout <- results[[3]][[2]][[2]]
corr_bias_y_holdout <- results[[4]][[2]][[2]]
holdout_list=list(good_y_holdout,corr_y_holdout,bias_y_holdout,corr_bias_y_holdout)

compare_list <- lapply(1:4, FUN=function(x){
  y_holdout <- holdout_list[[x]]
  fits <- fits_list[[x]]
  elpd_compare <- matrix(nrow=10,ncol=2)
  elpd_compare[1,] <- c(0,0)
  for (k in 2:10){
    if(k==2|k==4){
      bias = T
    } else {
      bias = F
    }
    elpd_compare[k,] <- elpd_comparer2(fit1=fits[[1]],fit2=fits[[j]],
                                       bias1=F,bias2=bias,y_holdout = y_holdout)
  }
  return(elpd_compare)
})

better_list <- lapply(compare_list,FUN=function(x){
  temp = rep("",10)
  temp[x[,1]+1.96*x[,2]<0] = "+"
  return(temp)
})
worse_list <- lapply(compare_list,FUN=function(x){
  temp = rep("",10)
  temp[x[,1]-1.96*x[,2]>0] = "-"
  return(temp)
})
better_list[[1]] <- rep("",10)
better_list[[2]] <- rep("",10)
worse_list[[1]] <- rep("",10)
worse_list[[2]] <- rep("",10)

palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
             "#CC79A7","#999999","#000000")

plotter2 <- function(elpds,title,M,eta,zeta,worse,better){
  elpds = as.data.frame(do.call(rbind,elpds))[c(1,4:9),]
  colnames(elpds) = c("Estimate","SE")
  to_plot <- mutate(elpds, Lower = Estimate - 1.96*SE, Upper = Estimate + 1.96*SE)
  Bayesbag_lab <- paste0("Opt. Bayesbag \n M = ",M)
  SMI_lab <- paste0("Opt. SMI \n \u03B7 = (",eta[1],", ",eta[2],")")
  coarse_lab <- paste0("Power Posterior \n (Coarsening) \n \u03B6 = ",zeta[1])
  labs <- c("Std. Bayes","Corr-Bias Bias","Bayesbag \n M = I",
            Bayesbag_lab,SMI_lab,coarse_lab,"Posterior Bootstrap")
  to_plot <- cbind(labs, data.frame(to_plot, row.names=NULL))
  colnames(to_plot) = c("Method",colnames(to_plot)[2:ncol(to_plot)])
  to_plot <- cbind(data.frame(to_plot, row.names=NULL),worse[c(1,4:9)])
  to_plot <- cbind(data.frame(to_plot, row.names=NULL),better[c(1,4:9)])
  colnames(to_plot) = c(colnames(to_plot)[1:(ncol(to_plot)-2)],"Worse","Better")
  ggplot(to_plot,aes(ymin = Lower, ymax = Upper, 
                     y=Estimate, x = fct_inorder(Method))) +
    theme_bw() + scale_colour_manual(values=palette) +
    geom_pointrange(size=1.1,fatten=4,aes(col=Method),show.legend = FALSE) +
    geom_text(aes(label = Worse),size = 6,
              nudge_y = 0.01, nudge_x = 0.1, color = "black") +
    geom_text(aes(label = Better),size = 6, 
              nudge_y = 0.01, nudge_x = 0.1, color = "black") +
    labs(title=element_blank(), y="MLPD",x=element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5))
}

plotnames <- list("Ex2 standard MLPDs","Ex2 correlated MLPDs",
                  "Ex2 biased MLPDs","Ex2 correlated-biased MLPDs")
ex2plots <- lapply(1:4,function(x){
  elpds <- elpds_list[[x]]
  M = results[[x]][[1]][[6]][[4]]
  eta = results[[x]][[1]][[7]][[3]]
  zeta = results[[x]][[1]][[8]][[3]]
  better = better_list[[x]]
  worse = worse_list[[x]]
  plot <- plotter2(elpds=elpds,M=M,eta=eta,zeta=zeta,better=better,worse=worse)
  ggsave(filename = paste0(plotnames[[x]],".pdf"),plot=plot,width=20,height=13,
         unit="cm",device=cairo_pdf)
})

