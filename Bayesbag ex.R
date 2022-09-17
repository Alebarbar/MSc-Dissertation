##Bayesbag ex 4
#This code will generate all the data & results (except graphs) for example 1
#N.B. this takes about 2 weeks to run on a 10-core intel 10900 and uses at least
#32GB of RAM, but is very parallel to should scale up to around 100 cores 
#It is set up to still be reproducible even if the different sections at the end
#"reproduce paper", "more reasonable model" etc. are run on different computers
#and the "skip" argument can be used to only calculate one of the three cases
#this allows the workload to be simply split over 12 machines if desired

#make sure you have rstan version 2.26 or later (this is later than the current 
#release version) as the stan syntax for arrays has changed and the code will 
#not run otherwise. Find instructions here:
#https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

#Setup
library(MASS)
library(rstan)
library(parallel)
library(rstudioapi)

#complex rng setup to allow reproducability even when running on 
#machines with different numbers of cores
RNGkind("L'Ecuyer-CMRG")
set.seed(38472693)
seedslist <- lapply(as.list(1:10), FUN = function(x){
  lapply(as.list(1:10), FUN = function(z){
    sample(1:10^9,1000,replace=T)
    })
  })

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

invlogit <- function(x){return(1/(1+exp(-x)))}

#function to generate synthetic data
HMELR <- function(I=3,J=8,K=100,D=3,reps=R,
                  beta=c(0.65,rep(1,D)),sigma_v=3,sigma_u=3,rho=0,bias=0){
  v = vector(length=K)
  u = matrix(nrow=J,ncol=K)
  d = c(sigma_v,rep(sigma_u,J))
  sh1 = diag(x=d^2,nrow=1+J,ncol=1+J)
  sh2 = matrix(rho,nrow=1+J,ncol=1+J)-diag(rep(rho,J+1))
  sh3 = matrix(d,nrow=1+J,ncol=1+J,byrow=T)
  sh4 = matrix(d,nrow=1+J,ncol=1+J,byrow=F)
  Sigma = sh1 + sh2*sh3*sh4
  for (k in 1:K){
    temp <- mvrnorm(mu=rep(0,1+J),Sigma=Sigma)
    v[k]=temp[1]
    u[,k]=temp[2:(J+1)]
    }
  Z = list()
  Z[[1]] = array(data=rnorm(n=I*J*K*D,mean=0,sd=1),dim=c(I,J,K,D))
  p = list()
  p[[1]] = array(dim=c(I,J,K))
  Y = array(dim=c(I,J,K))
  Y_holdout = list()
  for (i in 1:I){
    for (j in 1:J){
      for (k in 1:K){
        pre_p <- beta[1]+t(Z[[1]][i,j,k,])%*%beta[2:(D+1)]+u[j,k]+v[k]
        p[[1]][i,j,k] = invlogit(pre_p)
        Y[i,j,k] = rbinom(1,1,p[[1]][i,j,k]*(1-bias))
      }
    }
  }
  for (r in 1:reps){
    Y_holdout[[r]] = array(dim=c(I,J,K))
    Z[[1+r]] = array(data=rnorm(n=I*J*K*D,mean=0,sd=1),dim=c(I,J,K,D))
    p[[1+r]] = array(dim=c(I,J,K))
    for (i in 1:I){
      for (j in 1:J){
        for (k in 1:K){
          pre_p <- beta[1]+t(Z[[1+r]][i,j,k,])%*%beta[2:(D+1)]+u[j,k]+v[k]
          p[[1+r]][i,j,k] = invlogit(pre_p)
          Y_holdout[[r]][i,j,k] = rbinom(1,1,p[[1+r]][i,j,k]*(1-bias))
        }
      }
    }
  }
  params = list(v,u,Z,beta,p)
  return(list(Y,params,Y_holdout))
}

#function to use stan to sample from posteriors
stanner <- function(y, method="combo",rho=0, eta_v=1, eta_u=1, zeta=1,
                     warmup = warmup_iters, adapt_delta=0.99, max_treedepth=12,
                    chain_id = 1, chains = 10,iter=N, seed = sample(1:10^9,1),
                    pars = c("v","u","sigma_v","sigma_u","beta_0","beta")){
  stan_input=list(I,J,K,D,sigma_v_upper,
                  sigma_u_upper,beta_sd,rho,eta_v,eta_u,zeta,y)
  names(stan_input) <- c("I","J","K","D","sigma_v_upper","sigma_u_upper",
                         "beta_sd","rho","eta_v","eta_u","zeta","y")
  file = paste0("HMELR ", method, ".stan")
  model <- stan_model(file=file, model_name = method)
  fit <- sampling(model, data = stan_input, chains = chains, warmup = warmup,
                chain_id = chain_id, iter = iter, pars=pars, save_warmup=F,
                seed = seed, control = list("max_treedepth"=max_treedepth,
                                            "adapt_delta"=adapt_delta))
  return(fit)
}

#function to calculate elpd of stan fitted model using holdout data
mcer <- function(cnst,total_sd,draws=100){
  temp <- sapply(1:draws,FUN=function(x){
    Z <- rnorm(1,1)
    p = invlogit(Z*total_sd+cnst)
    q = invlogit(-Z*total_sd-cnst)
    return(c(p,q))
  })
  p = mean(temp[1,])
  q = mean(temp[2,])
  return(c(p,q))
}

elpdhelper <- function(input,reps,I,J,K,y_holdout,mc_draws,bias=0,supplied_p=NA){
  log_pd <- vector(length=reps*I*J*K)
  beta_0 = input[[3]]
  total_sd = input[[4]]
  acc = 1
  for (r in 1:reps){
    for (k in 1:K){
      v = input[[1]][k]
      for (j in 1:J){
        u = input[[2]][j,k]
        cnst = beta_0 + v + u
        if (anyNA(supplied_p)){ #generate p unless supplied
          temp = mcer(cnst=cnst,total_sd=total_sd, draws=mc_draws)
          p = temp[1]*(1-bias)
          q = temp[2]*(1-bias)+bias
        }
        for (i in 1:I){
          if (!anyNA(supplied_p)){
            p = supplied_p[[r+1]][i,j,k]*(1-bias)
            q = 1-p
          }
          if (y_holdout[[r]][i,j,k] == 1){
            if (p == 1){log_pd[acc] <- log1p(-q)}
            else{log_pd[acc] <- log(p)}
          }
          else{
            if (q == 1){log_pd[acc] <- log1p(-p)}
            else{log_pd[acc] <- log(q)}}
          acc = acc + 1
        }
      }
    }
  }
  return(log_pd)
}

elpder <- function(fit, y_holdout,mc_draws=100,cl=NA,return_lpds=F){
  draws <- rstan::extract(fit, pars = c("v","u","beta_0","beta"))
  I = length(y_holdout[[1]][,1,1])
  J = length(y_holdout[[1]][1,,1])
  K = length(y_holdout[[1]][1,1,])
  S = dim(draws[[1]])[1]
  reps <- length(y_holdout)
  
  input <- lapply(1:S,FUN=function(s){
    v <- draws[[1]][s,]
    u <- draws[[2]][s,,]
    beta_0 <- draws[[3]][s]
    beta <- draws[[4]][s,]
    total_sd = sqrt(sum(beta^2))
    return(list(v,u,beta_0,total_sd))
  })
  if(anyNA(cl)){
    log_pd = lapply(X=input,FUN = elpdhelper,reps=reps,I=I,J=J,K=K,
                    y_holdout=y_holdout, mc_draws=mc_draws)
  }else{
    log_pd = parLapply(cl,X=input,fun = elpdhelper,reps=reps,I=I,J=J,K=K,
                       y_holdout=y_holdout, mc_draws=mc_draws)
  }
  log_pd = do.call(rbind,log_pd)
  lpds <- matrixStats::colLogSumExps(log_pd) - log(nrow(log_pd))
  elpd <- mean(lpds)
  se <- sd(lpds)/sqrt(length(lpds))
  if(return_lpds){
    return(list(c(elpd,se),lpds))
  }else{
    return(c(elpd,se))
  }
}

#function to extract information about the posterior 
#parameter draws from the fitted stan model
post_paramer <- function(fit){
  post_params_summary <- summary(fit)$summary
  post_param_quantiles = matrix(nrow=0,ncol=fidelity)
  acc = 1
  for (x in rstan::extract(fit)){
    add <- prod(dim(x)[-1])
    temp = split(x,rep(acc:(acc+add-1), each = dim(x)[1]))
    temp = lapply(temp,quantile,probs=(0:(fidelity-1))/(fidelity-1))
    temp = lapply(temp,matrix,ncol=fidelity)
    temp = do.call(rbind,temp)
    post_param_quantiles <- rbind(post_param_quantiles,temp)
    acc = acc + add
  }
  percentages = paste0(round(100*(0:(fidelity-1))/(fidelity-1),digits=2),"%")
  colnames(post_param_quantiles) = percentages
  rownames(post_param_quantiles) = rownames(post_params_summary)
  return(list(post_params_summary,post_param_quantiles))
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

#function to calculate the posteriors and resulting elpds & parameter estimates
#for each technique for a given dataset
suite <- function(y_gen,y_gen_holdout,seeds = sample(1:10^9,1000)){
  #setup cluster to use throughout
  cl <- makeCluster(cluster_size)
  clusterExport(cl,varlist = list("y_gen","y_gen_holdout","seeds"),
                envir=environment())
  clusterExport(cl,varlist = list("sigma_v_upper","sigma_u_upper","beta_sd","I",
                "J","K","D","N","chain_num","R","warmup_iters","N_Bag",
                "Bag_chains","grid_scale","coarse_scale","fidelity","invlogit",
                "stanner","mcer","elpdhelper","elpder","post_paramer",
                "diagnoser"))
  clusterEvalQ(cl, {
    library(rstan)
    library(rstudioapi)
    library(loo)
    options(mc.cores = 10)
  })
  clusterSetRNGStream(cl,iseed=seeds[10])
  
  print("running std Bayes and fixed-rho")
  
  index = list()
  #standard Bayesian posterior, no rho model 
  index[[1]] = list("std Bayes",0)
  #Bayesian posterior with fixed rho=0.99
  index[[2]] = list("fixed-rho Bayes",0.99)
  
  holder <- parLapply(cl,index, fun=function(x){ 
    fit <- stanner(y=y_gen,rho=x[[2]], chains = chain_num)
    post_params <- post_paramer(fit)
    elpd = elpder(fit, y_holdout = y_gen_holdout)
    diagnostics <- diagnoser(fit)
    remove(fit)
    gc()
    return(list(elpd,post_params,diagnostics))
  })
  
  std_bayes_elpd <- holder[[1]][[1]]
  std_bayes_post_params <- holder[[1]][[2]]
  std_bayes_diagnostics <- holder[[1]][[3]]
  fixed_rho_elpd <- holder[[2]][[1]]
  fixed_rho_post_params <- holder[[2]][[2]]
  fixed_rho_diagnostics <- holder[[2]][[3]]
  remove(holder)
  print("finished std and fixed rho")
  gc()
  
  #Bayesbag model with B samples
  
  if(B==0){ #skip Bayesbag if B=0
    print("skipping Bayesbag")
    bayesbag_post_params <- NA
    bayesbag_elpd <- NA
    bayesbag_diagnostics <- NA
  }else{
  
  print("starting Bayesbag")
  clusterEvalQ(cl, {
     options(mc.cores = Bag_chains)
   })
  bayesbag_sflist <- parLapplyLB(cl,(1:B)*Bag_chains,fun = function(cid){
    set.seed(seeds[100+cid/Bag_chains])
    gc()
    y_bagged <- array(dim=c(I,J,K))
    for (j in 1:J){
      for (k in 1:K){
        y_bagged[,j,k] = sample(x=y_gen[,j,k],size=I,replace=T)
      }
    }
    return(stanner(y=y_bagged,iter=N_Bag,chain_id=cid, 
                   warmup = N_Bag/2, chains = Bag_chains, 
                   seed = seeds[100+cid/Bag_chains]))
  })
  bagged_diagnostics <- lapply(bayesbag_sflist, diagnoser)
  bayesbag_fit <- sflist2stanfit(bayesbag_sflist)
  remove(bayesbag_sflist)
  bayesbag_post_params <- post_paramer(bayesbag_fit)
  bayesbag_diagnostics <- diagnoser(bayesbag_fit)
  clusterSetRNGStream(cl,iseed=seeds[301])
  bayesbag_elpd <- elpder(bayesbag_fit, y_holdout = y_gen_holdout,cl=cl)
  remove(bayesbag_fit)
  gc()
  print("Bayesbag done")
  }
  
  #grid-search SMI
  
  if(grid_scale==0){ #skip SMI if grid_scale=0
    print("skipping grid_search SMI")
    SMI_elpd_vals <- NA
    SMI_opt_elpd <- NA
    SMI_opt_post_params <- NA
    SMI_opt_diagnostics <- NA
    eta_opt <- NA
  }else{
  
  print("starting grid-search SMI")
  grid <- list()
  for (p in 0:(grid_scale-1)){
    for (q in 0:(grid_scale-1)){
    grid[[p*grid_scale+q+1]] = c(p,q)
    }
  }
  holder <- parLapplyLB(cl,grid,fun=function(coords){
    gc()
    set.seed(seeds[500+p*grid_scale+q+1])
    p = coords[1]
    q = coords[2]
    eta_v=p/(grid_scale-1)
    eta_u=q/(grid_scale-1)
    SMI_fit <- stanner(y=y_gen,iter=N_Bag, warmup=N_Bag/2,
                       eta_v=eta_v,eta_u=eta_u,
                       chain_id = (p*grid_scale+q+1)*Bag_chains, 
                       chains=Bag_chains, seed = seeds[700+p*grid_scale+q+1])
    SMI_post_params <- post_paramer(SMI_fit)
    SMI_diagnostics <- diagnoser(SMI_fit)
    temp <- elpder(SMI_fit,y_holdout = y_gen_holdout,return_lpds =T)
    elpd <- temp[[1]]
    lpds <- temp[[2]]
    remove(SMI_fit)
    return(list(elpd,SMI_post_params,SMI_diagnostics,lpds))
  })
  
  SMI_elpd_vals <- array(dim=c(grid_scale,grid_scale,2))
  for (p in 0:(grid_scale-1)){
    for (q in 0:(grid_scale-1)){
      SMI_elpd_vals[p+1,q+1,] = holder[[p*grid_scale+q+1]][[1]]
    }
  }
  
  SMI_opt_coords = which(SMI_elpd_vals[,,1]==max(SMI_elpd_vals[,,1]),arr.ind=T)
  trans_coords = (SMI_opt_coords[1]-1)*grid_scale+SMI_opt_coords[2]
  t1 <- holder[[trans_coords]][[4]]
  t2 <- holder[[grid_scale^2]][[4]]
  t3 = t1 - t2
  t4 = mean(t3) - 1.96*sd(t3)/sqrt(length(t3))
  if(t4<=0){
    SMI_opt_coords = c(grid_scale,grid_scale)
    trans_coords = (SMI_opt_coords[1]-1)*grid_scale+SMI_opt_coords[2]
  }
  SMI_opt_elpd <- holder[[trans_coords]][[1]]
  SMI_opt_post_params <- holder[[trans_coords]][[2]]
  SMI_opt_diagnostics <- holder[[trans_coords]][[3]]
  eta_opt <- (SMI_opt_coords-1)/(grid_scale-1)
  remove(holder,t1,t2,t3,t4)
  print("finished grid-serach SMI")
  }
  
  #Coarsening
  if(coarse_scale==0){ #skip Coarsening if coarse_scale=0
    print("skipping Coarsening")
    coarse_elpds <- NA
    coarse_opt_elpd <- NA
    coarse_opt_post_params <- NA
    coarse_opt_diagnostics <- NA
    zeta_opt <- NA
  }else{
    print("starting Coarsening")
    holder <- parLapplyLB(cl,(0:(coarse_scale-1)),fun=function(x){
      set.seed(seeds[900+x])
      zeta <- x/(coarse_scale-1)
      fit <- stanner(y=y_gen, iter=N_Bag, warmup=N_Bag/2,
                         zeta=zeta, chain_id = x*Bag_chains, 
                         chains=Bag_chains, seed = seeds[950+x])
      post_params <- post_paramer(fit)
      diagnostics <- diagnoser(fit)
      temp <- elpder(fit,y_holdout = y_gen_holdout,return_lpds = T)
      elpds <- temp[[1]]
      lpds <- temp[[2]]
      remove(fit)
      return(list(elpds,post_params,diagnostics,lpds))
    })
  
    coarse_elpds <- matrix(ncol=2,nrow=coarse_scale)
    for (x in 1:coarse_scale){
      coarse_elpds[x,] = holder[[x]][[1]]
    }
    coarse_opt_coord = which(coarse_elpds[,1]==max(coarse_elpds[,1]),arr.ind=T)
    t1 <- holder[[coarse_opt_coord]][[4]]
    t2 <- holder[[coarse_scale]][[4]]
    t3 <- t1 - t2
    t4 <- mean(t3) - 1.96*sd(t3)/sqrt(length(t3))
    if(t4 <= 0){
      coarse_opt_coord <- coarse_scale
    }
    coarse_opt_elpd <- holder[[coarse_opt_coord]][[1]]
    coarse_opt_post_params <- holder[[coarse_opt_coord]][[2]]
    coarse_opt_diagnostics <- holder[[coarse_opt_coord]][[3]]
    zeta_opt <- (coarse_opt_coord-1)/(coarse_scale-1)
    remove(holder,t1,t2,t3,t4)
    print("finished Coarsening")
  }
    
  stopCluster(cl)
  
  elpd_list <-list(std_bayes_elpd, fixed_rho_elpd, bayesbag_elpd, SMI_opt_elpd,
                  eta_opt,SMI_elpd_vals, coarse_opt_elpd, zeta_opt,coarse_elpds)
  post_param_list <-list(std_bayes_post_params,
                         fixed_rho_post_params,bayesbag_post_params,
                         SMI_opt_post_params,coarse_opt_post_params)
  diagnostic_list <- list(std_bayes_diagnostics, fixed_rho_diagnostics,
                          bayesbag_diagnostics, bagged_diagnostics, 
                          SMI_opt_diagnostics, coarse_opt_diagnostics)
  return(list(elpd_list,post_param_list,diagnostic_list))
}

#function to generate data and apply all the tests
totaliser <- function(data_good=NA,data_miss=NA,
                      data_bias=NA, skip = c(F,F,F),sigma_v=3,sigma_u=3){
  namelist1 <- c("std Bayes","fixed-rho Bayes","Bayesbag","Optimal SMI",
                 "Optimal Coarsening")
  namelist2<- c("std Bayes","fixed-rho Bayes","Bayesbag","Optimal SMI",
                "Opt SMI Etas","all SMI elpds","Optimal Coarsening",
                "Opt Coarsening zeta","all Coarsening elpds")
  namelist3 <- c("std Bayes","fixed-rho Bayes","Bayesbag","Bagged",
                 "Optimal SMI","Optimal Coarsening")
  elpd_NAs <- rep(NA,9)
  post_param_NAs <- rep(NA,5)
  diagnostics_NAs <- rep(NA,6)
                
  set.seed(seeds_to_use[[1]][1])
  if(anyNA(data_good)){
    data_good = HMELR(I=I,J=J,K=K,D=D,
                      beta = c(0.65,rep(1,D)),sigma_v=sigma_v,sigma_u=sigma_u)
  }
  y_good = data_good[[1]]
  params_good = data_good[[2]]
  y_good_holdout = data_good[[3]]
  if(skip[1]){
    good_elpds = elpd_NAs
    good_post_params = post_param_NAs
    good_diagnostics = diagnostics_NAs
  }else{
    print("starting good data eval") 
    output <- suite(y_gen=y_good,
                    y_gen_holdout = y_good_holdout,seeds=seeds_to_use[[1]])
    good_elpds <- output[[1]]
    good_post_params <- output[[2]]
    good_diagnostics <- output[[3]]
    names(good_elpds) = namelist2
    names(good_post_params) = namelist1
    good_diagnostics <- output[[3]]
    names(good_diagnostics) = namelist3
    print("good data eval done")
    gc()
  }
  
  set.seed(seeds_to_use[[2]][1])
  if(anyNA(data_miss)){
    data_miss = HMELR(I=I,J=J,K=K,D=D,beta = c(0.65,rep(1,D)),
                      rho = 0.99,sigma_v=sigma_v,sigma_u=sigma_u)
  }
  y_miss = data_miss[[1]]
  params_miss = data_miss[[2]]
  y_miss_holdout = data_miss[[3]]
  if(skip[2]){
    miss_elpds = elpd_NAs
    miss_post_params = post_param_NAs
    miss_diagnostics = diagnostics_NAs
  }else{
    print("starting miss data eval") 
    output <- suite(y_gen=y_miss,
                    y_gen_holdout = y_miss_holdout,seeds=seeds_to_use[[2]])
    miss_elpds <- output[[1]]
    miss_post_params <- output[[2]]
    names(miss_elpds) = namelist2
    names(miss_post_params) = namelist1
    miss_diagnostics <- output[[3]]
    names(miss_diagnostics) = namelist3
    print("miss data eval done")
    gc()
  }
  
  set.seed(seeds_to_use[[3]][1])
  if(anyNA(data_bias)){
    data_bias = HMELR(I=I,J=J,K=K,D=D,beta = c(0.65,rep(1,D)),
                      bias = 0.1,sigma_v=sigma_v,sigma_u=sigma_u)
  }
  y_bias = data_bias[[1]]
  params_bias = data_bias[[2]]
  y_bias_holdout = data_bias[[3]]
  if(skip[3]){
    bias_elpds = elpd_NAs
    bias_post_params = post_param_NAs
    bias_diagnostics = diagnostics_NAs
  }else{
    print("starting bias data eval") 
    output <- suite(y_gen=y_bias,
                    y_gen_holdout = y_bias_holdout,seeds=seeds_to_use[[3]])
    bias_elpds <- output[[1]]
    bias_post_params <- output[[2]]
    names(bias_elpds) = namelist2
    names(bias_post_params) = namelist1
    bias_diagnostics <- output[[3]]
    names(bias_diagnostics) = namelist3
    print("bias data eval done")
    gc()
  }
  
  elpds <- list(good_elpds,miss_elpds,bias_elpds)
  data <- list(data_good,data_miss,data_bias)
  post_params <- list(good_post_params,miss_post_params,bias_post_params)
  diagnostics <- list(good_diagnostics, miss_diagnostics, bias_diagnostics)
  return(list(elpds,data,post_params,diagnostics))
}

fidelity=51 #setting for how much detail to record for post param dists
Bag_chains=2
N_Bag = 2000
N=2000
chain_num = 4
results <- list()

#Section 1: Replicate paper
seeds_to_use <- seedslist[[1]]
sigma_v_upper = 100
sigma_u_upper = 100
beta_sd = 10
I=3
J=8
K=100
D=3
R = 4
warmup_iters=N/2
cluster_size = parallel::detectCores()/Bag_chains
B=100
grid_scale = 11
coarse_scale = 11
std_data = totaliser(skip = c(T,T,T))[[2]]

#after this point you can only include the section you want to run
#And use the 'skip' argument to skip running some of the data sections

results[[1]] <- totaliser(skip = c(F,F,F), data_good = std_data[[1]],
                          data_miss = std_data[[2]],data_bias = std_data[[3]])
to_save <- results[[1]]
save(list = "to_save",file = "replication_results.RData")
remove(to_save)


#Section 2: Concentrated priors
#try less diffuse priors
N = 5000
chain_num = 10
seeds_to_use <- seedslist[[2]]
sigma_v_upper = 10
sigma_u_upper = 10
beta_sd = 1
I=3
J=8
K=100
D=3
R = 4
warmup_iters=N/2
cluster_size = parallel::detectCores()/Bag_chains
B=100
grid_scale = 11
coarse_scale = 11

results[[2]] <- totaliser(skip = c(F,F,F),data_good = std_data[[1]],
                          data_miss = std_data[[2]],data_bias = std_data[[3]])
to_save <- results[[2]]
save(list = "to_save",file = "conc_results.RData")
remove(to_save)

#Section 3: More data model
#See if more data can also overpower prior
N = 2000
chain_num = 4
seeds_to_use <- seedslist[[3]]
sigma_v_upper = 100
sigma_u_upper = 100
beta_sd = 10
I=30
J=30
K=30
D=3
R = 2
warmup_iters=N/2
cluster_size = parallel::detectCores() #not dividing by bag_chains speeds up 
#elpd eval so use this when not running Bayesbag, SMI, or Coarsening
B=0 #skipped here
grid_scale <- 0 #skipped here
coarse_scale = 0 #skipped here

results[[3]] <- totaliser(skip = c(F,F,T))
to_save <- results[[3]]
save(list = "to_save",file = "more_data_results.RData")
remove(to_save)

save(list="results",file="full_results.RData")

