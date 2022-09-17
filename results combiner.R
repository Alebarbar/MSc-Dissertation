#Takes any combination of individual results files, accepting each 
#test either as one file or split into three for the good, miss and bias results
#(If split into three it will first combine them into a single file and save this)
#loads 'results' into memory and saves a 'full_results.RData' file

loader <- function(name){
  NAlist <- list(rep(list(rep(NA,5)),3),list(NA,NA,NA),rep(list(rep(NA,4)),3))
  fname <- paste0(name,"_results.RData")
  if(!file.exists(fname)){
    goodname <- paste0(name,"_results_good.RData")
    if(file.exists(goodname)){
      load(goodname)
      good <- to_save
      remove(to_save)
    } else{
      good = NAlist
    }
    
    missname <- paste0(name,"_results_miss.RData")
    if(file.exists(missname)){
      load(missname)
      miss <- to_save
      remove(to_save)
    }else{
      miss = NAlist
    }
    
    biasname <- paste0(name,"_results_bias.RData")
    if(file.exists(biasname)){
      load(biasname)
      bias <- to_save
      remove(to_save)
    }else{
      bias = NAlist
    }
    
    elpds <- list(good[[1]][[1]],miss[[1]][[2]],bias[[1]][[3]])
    names(elpds) <- c("good_elpds","miss_elpds","bias_elpds")
    
    if(!anyNA(good[[2]])){
      data <- good[[2]]
    }else if(!anyNA(miss[[2]])){
      data <- miss[[2]]
    }else{
      data <- bias[[2]]
    }
    names(data) <- c("data_good","data_miss","data_bias")
    
    post_params <- list(good[[3]][[1]],miss[[3]][[2]],bias[[3]][[3]])
    names(post_params) <- c("good_post_params","miss_post_params","bias_post_params")
    
    to_save <- list(elpds,data,post_params)
    save(to_save,file=fname)
    remove(to_save)
  }
  load(fname)
  return(to_save)
}

if(file.exists("full_results.RData")){  
  load("full_results.RData")
  }else{
  results <- lapply(c("replication","conc","more_data"),loader)
  save(list = "results", file = "full_results.RData")
  }
