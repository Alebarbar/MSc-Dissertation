library(ggplot2)
library(tidyverse)
library(cowplot)
library(Cairo)

invlogit <- function(x){return(1/(1+exp(-x)))}

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


palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
             "#0072B2", "#D55E00", "#CC79A7","#999999")

plotter <- function(d,true_elpd,title,eta,zeta){
  to_plot <- as.data.frame(t(as.data.frame(d[c(1:4,7)])))
  colnames(to_plot) = c("Estimate","SE")
  to_plot <- mutate(to_plot, Lower = Estimate - 1.96*SE, Upper = Estimate + 1.96*SE)
  SMI_lab <- paste0("Opt. SMI \n \u03B7 = (",eta[1],", ",eta[2],")")
  coarse_lab <- paste0("Opt. Power Posterior \n (Coarsening) \n \u03B6 = ",zeta[1])
  labs <- c("Std. Bayes","Fixed-Rho Bayes","BayesBag \n M = I",SMI_lab,coarse_lab)
  to_plot <- cbind(labs, data.frame(to_plot, row.names=NULL))
  colnames(to_plot) = c("Method",colnames(to_plot)[2:5])
  plot <- ggplot(to_plot,aes(ymin=Lower, ymax=Upper,
                             y=Estimate, x=fct_inorder(Method))) +
    theme_bw() + scale_colour_manual(values=palette) +
    geom_pointrange(size=1.1,fatten=4,aes(col=Method),show.legend = FALSE) +
    labs(title=title, y="MLPD",x=element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5))
  scale <- layer_scales(plot)$y$range$range
  s <- (scale[2]-scale[1])*0.03
  y_min <- scale[1]-(scale[2]-scale[1])*0.25
  if(log(1/2)>y_min){
    plot = plot + geom_hline(yintercept=log(1/2), color="red") + 
    annotate("text", x=0.7, y=log(1/2)-s, label="-0.69",colour="red")
  }
  y_max <- scale[2]+(scale[2]-scale[1])*0.25
  if(true_elpd<y_max){
    plot = plot +geom_hline(yintercept=true_elpd, color="blue") + 
    annotate("text", x=0.7, y=true_elpd-s, label=round(true_elpd,2),colour="blue")
  }else{
    plot = plot + annotate("text", x=0.7, y=y_max-s, 
                           label=paste0("^ ",round(true_elpd,2)),colour="blue")
  }
  return(plot)
}

visualiser <- function(d){
  temp = t(rbind(d,(0:50)/50))
  colnames(temp) = c("val","per")
  temp = as.data.frame(temp)
  temp = mutate(temp,den=100/(lead(val,default=1000)-val))
  ggplot(as.data.frame(temp),aes(y=den,x=val))+geom_point()
}

#load("full_results.RData")

rep_true_elpds <- lapply(1:3,FUN = function(x){
  if(x==3){
    bias = 0.1
  }else{
    bias=0
  }
  y_holdout = results[[1]][[2]][[x]][[3]]
  R = length(y_holdout)
  dims = dim(y_holdout[[1]])
  I = dims[1]
  J = dims[2]
  K = dims[3]
  data <- results[[1]][[2]]
  v=data[[x]][[2]][[1]]
  u=data[[x]][[2]][[2]]
  beta_0=data[[x]][[2]][[4]][1]
  total_sd=sqrt(sum(data[[x]][[2]][[4]][-1]^2))
  
  p_list = data[[x]][[2]][[5]]
  int_input <- list("v"=v,"u"=u,"beta_0"=beta_0,"total_sd"=total_sd)
  inted <- mean(elpdhelper(input=int_input,reps=R,I=I,J=J,K=K,
                         y_holdout=y_holdout,bias=bias, mc_draws=1000))
  noint <- mean(elpdhelper(input=int_input,reps=R,I=I,J=J,K=K,
                y_holdout=y_holdout, mc_draws=1000,bias=bias, supplied_p = p_list))
  return(c(inted,noint))
})

conc_true_elpds <- lapply(1:3,FUN = function(x){
  if(x==3){
    bias = 0.1
  }else{
    bias=0
  }
  y_holdout = results[[2]][[2]][[x]][[3]]
  R = length(y_holdout)
  dims = dim(y_holdout[[1]])
  I = dims[1]
  J = dims[2]
  K = dims[3]
  data <- results[[2]][[2]]
  v=data[[x]][[2]][[1]]
  u=data[[x]][[2]][[2]]
  beta_0=data[[x]][[2]][[4]][1]
  total_sd=sqrt(sum(data[[x]][[2]][[4]][-1]^2))
  p_list = data[[x]][[2]][[5]]
  int_input <- list("v"=v,"u"=u,"beta_0"=beta_0,"total_sd"=total_sd)
  inted <- mean(elpdhelper(input=int_input,reps=R,I=I,J=J,K=K,
                           y_holdout=y_holdout,bias=bias, mc_draws=1000))
  noint <- mean(elpdhelper(input=int_input,reps=R,I=I,J=J,K=K,
                           y_holdout=y_holdout, mc_draws=1000,
                           bias=bias, supplied_p = p_list))
  return(c(inted,noint))
})

rep_good_elpds <- results[[1]][[1]][[1]]
rep_miss_elpds <- results[[1]][[1]][[2]]
rep_bias_elpds <- results[[1]][[1]][[3]]

conc_good_elpds <- results[[2]][[1]][[1]]
conc_miss_elpds <- results[[2]][[1]][[2]]
conc_bias_elpds <- results[[2]][[1]][[3]]

more_data_good_elpds <- results[[3]][[1]][[1]]
more_data_miss_elpds <- results[[3]][[1]][[2]]
more_data_bias_elpds <- results[[3]][[1]][[3]]

#generate replication mlpds plot
rep_good_plot <- plotter(rep_good_elpds,title="Uncorrelated Case",
             true_elpd = rep_true_elpds[[1]][1],eta=rep_good_elpds[[5]],
             zeta = rep_good_elpds[[8]])
rep_miss_plot <- plotter(rep_miss_elpds,title="Correlated Case",
                  true_elpd = rep_true_elpds[[2]][1],eta=rep_miss_elpds[[5]],
                         zeta = rep_miss_elpds[[8]])
rep_plot_1 <- plot_grid(rep_good_plot,rep_miss_plot)
ggsave("rep_plot_1.pdf",rep_plot_1,width=30,height=15,unit="cm",device=cairo_pdf)

rep_plot_2 <- plotter(rep_bias_elpds,title="Biased Case",
           true_elpd = rep_true_elpds[[3]][1],eta=rep_bias_elpds[[5]],
                      zeta = rep_bias_elpds[[8]])
ggsave("rep_plot_2.pdf",rep_plot_2,width=15,height=15,unit="cm",device=cairo_pdf)

#generate conc mlpds plot
conc_good_plot <- plotter(conc_good_elpds,title="Uncorrelated Case",
                         true_elpd = conc_true_elpds[[1]][1],eta=conc_good_elpds[[5]],
                         zeta = conc_good_elpds[[8]])
conc_miss_plot <- plotter(conc_miss_elpds,title="Correlated Case",
                         true_elpd = conc_true_elpds[[2]][1],eta=conc_miss_elpds[[5]],
                         zeta = conc_miss_elpds[[8]])
conc_plot_1 <- plot_grid(conc_good_plot,conc_miss_plot)
ggsave("conc_plot_1.pdf",conc_plot_1,width=30,height=15,unit="cm",device=cairo_pdf)

conc_plot_2 <- plotter(conc_bias_elpds,title="Biased Case",
                      true_elpd = conc_true_elpds[[3]][1],eta=conc_bias_elpds[[5]],
                      zeta = conc_bias_elpds[[8]])
ggsave("conc_plot_2.pdf",conc_plot_2,width=15,height=15,unit="cm",device=cairo_pdf)