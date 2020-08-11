# paths
fig_path = "/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/figures"
options("scipen"=100, "digits"=4)

# libs
library(tidyverse)
library(cowplot)
library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(readxl)
library(shinystan)

# ggplot2 theme
theme_custom = theme_light() + theme(text = element_text(size = 8),
                                     axis.text=element_text(size=7),
                                     strip.text=element_text(size=8))
theme_set(theme_custom)

# functions
qsum = function(x) c(mean=mean(x),quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975)))
qsum2 = function(x) c(mean=mean(x),quantile(x,probs=c(0.025,0.5,0.975)))
distsum = function(x,dist) {
  if(dist=="gamma") r = paste0(round(x[1]/x[2],2)," (",round(sqrt(x[1]/(x[2]^2)),2),")")
  if(dist=="beta") r = paste0(round(x[1]/(x[1]+x[2]),2)," (",round(sqrt(x[1]*x[2]/( (x[1] + x[2])^2*(x[1]+x[2]+1))),2),")")
  return(r)
}
get_gammapar = function(mu,sigma) {
  alpha = mu^2/sigma^2
  beta = mu/sigma^2
  return(c(alpha=alpha,beta=beta))
}
get_betapar = function(mu,sigma) {
  alpha = -(mu^3 + mu*sigma^2 - mu^2)/sigma^2
  beta = (mu^3 + (sigma^2 + 1)*mu - 2*mu^2 - sigma^2)/sigma^2
  return(c(alpha=alpha,beta=beta))
}
get_lnormpar = function(mu,sigma) {
  logmu = log(mu/sqrt(1+sigma^2/mu^2)); 
  logsigma = sqrt(log(1+sigma^2/mu^2));
  return(c(logmu=logmu,logsigma=logsigma))
}

check_gamma = function(shape,rate) {
  xx = rgamma(1000,shape,rate)
  plot(density(xx),main=paste0("~Gamma(",shape,",",rate,")"))
  cat(paste0("Distribution of rate: ",round(mean(xx),2)," [",round(quantile(xx,0.025),2),"-",round(quantile(xx,0.975),2),"] years^-1\n"))
  cat(paste0("Distribution of duration: ",round(mean(1/xx),2)," [",round(quantile(1/xx,0.025),2),"-",round(quantile(1/xx,0.975),2),"] years"))
}
lognormal_repar = function(raw,p) {
  Lp=c(NA,NA)
  Lp[1] = log(p[1]/sqrt(1+p[2]^2/p[1]^2))
  Lp[2] = sqrt(log(1+p[2]^2/p[1]^2))
  repar = exp(Lp[1] + raw * Lp[2])
  return(repar)
}
lognormal_getpar = function(p) {
  Lp=c(NA,NA)
  Lp[1] = log(p[1]/sqrt(1+p[2]^2/p[1]^2))
  Lp[2] = sqrt(log(1+p[2]^2/p[1]^2))
  return(Lp)
}
check_lognorm = function(m,s) {
  xx = rnorm(1000,0,1)
  xx = lognormal_repar(xx,c(m,s))
  plot(density(xx),main=paste0("~Lognormal with mean ",m," and sd ",s,""))
  cat(paste0("Distribution of rate: ",round(mean(xx),2)," [",round(quantile(xx,0.025),2),"-",round(quantile(xx,0.975),2),"] years^-1\n"))
  cat(paste0("Distribution of duration: ",round(mean(1/xx),2)," [",round(quantile(1/xx,0.025),2),"-",round(quantile(1/xx,0.975),2),"] years"))
}
cortocov = function(corr,sd) sd %*% t(sd) * corr

logit = function(x) log(x/(1-x))
inv.logit = function(x) exp(x)/(1+exp(x))

g_legend = function(a.gplot) { 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 
log_growth = function(t,tau,delay,slope) tau / (1+exp(-slope*(t+1-delay)))

# generate data list for stan
generate_base_D = function(indic_data,survey_data,key,starting_year,S,C,D,E,G,R=2) {
  tmp_indic_data = filter(indic_data,year>starting_year,year<=2018,Country.Name2 %in% key$Country.Name2) %>% arrange(CountryID) %>% mutate(year2=year-starting_year)
  tmp_survey_data = filter(survey_data,year>starting_year,year<=2018,Country.Name2 %in% key$Country.Name2) %>% arrange(CountryID) %>% mutate(year2=year-starting_year)
  if(G==1) tmp_survey_data$CountryID = 1
  tmp_suppression_data = filter(tmp_indic_data,!is.na(control_prop))
  if(G==1) tmp_suppression_data$CountryID = 1
  D = list(
    starting_year=starting_year,
    # controls
    inference=0,
    # simulation
    t0=0,
    S=S,
    ts=1:S,
    C=C,
    D=D,
    E=E,
    G=G,
    R=R,
    L=length(unique(tmp_indic_data$year)),
    M=dim(tmp_survey_data)[[1]],
    N=dim(tmp_suppression_data)[[1]],
    # initial values
    init_pop=filter(indic_data,year==starting_year,Country.Name %in% key$Country.Name) %>% arrange(CountryID) %>% pull(adult_pop) %>% structure(.,.Dim=G),
    init_prev=filter(indic_data,year==starting_year,Country.Name %in% key$Country.Name) %>% arrange(CountryID) %>% pull(adult_prevalence) %>% structure(.,.Dim=G),
    rollout=as.array(key$Rollout-starting_year),
    life_expectancy=key$life_expectancy,
    # indicator data (A-E)
    A_est=pull(tmp_indic_data,adult_prevalence) %>% matrix(.,nrow=G,byrow=TRUE),
    B_est=pull(tmp_indic_data,adult_art_coverage) %>% matrix(.,nrow=G,byrow=TRUE),
    C_est=pull(tmp_indic_data,adult_mortality) %>% matrix(.,nrow=G,byrow=TRUE),
    D_est=pull(tmp_indic_data,adult_pop) %>% matrix(.,nrow=G,byrow=TRUE),
    E_est=pull(tmp_indic_data,pmtct) %>% matrix(.,nrow=G,byrow=TRUE),
    # survey data (F)
    F_country = tmp_survey_data$CountryID,
    F_region = tmp_survey_data$Reg,
    F_n = tmp_survey_data$n_geno,
    F_k = tmp_survey_data$n_nnrti,
    F_t = tmp_survey_data$year-starting_year,
    # suppression data (G)
    G_country=pull(tmp_suppression_data,CountryID),
    G_est=pull(tmp_suppression_data,control_prop),
    G_t=pull(tmp_suppression_data,year2)
  )
  return(D)
}




# generate bash files for CMDSTAN
bashfile_rdump = function(model_name,data_file,warmup=500,iter=100,adapt_delta=0.8,max_depth=12,init=0,type="cmdstan",timelimit=96,chains=10) {
  if(type=="cmdstan") {
    data_file$inference=0
    with(data_file,stan_rdump(ls(data_file),file=paste0("models/data_SIM_",model_name,".R")))
    data_file$inference=1
    with(data_file,stan_rdump(ls(data_file),file=paste0("models/data_S_",model_name,".R")))
    if(length(init)!=1) {
      with(init,stan_rdump(ls(init),file=paste0("models/init_",model_name,".R")))
      init = paste0("init_",model_name,".R")
    }
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=",data_file$G),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_SIM_",model_name,".R output file=SIM_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10"),
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,".R output file=S_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
                  )
    writeLines(tt_launch,paste0("models/sb_",model_name,".sh"))
  }
  if(type=="rstan") {
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  "module load R",
                  paste0("Rscript run_",model_name,".R"))
    writeLines(tt_launch,paste0("sbr_",model_name,".sh"))
    tt_Rfile = c("library(rstan);options(mc.cores = parallel::detectCores())",
                 paste0("load('data_",model_name,".Rdata')"),
                 paste0("D_",model_name,"$inference = 0"),
                 paste0("SIM_",model_name," = stan('",model_name,".stan',data=D_",model_name,",chains=",chains,",iter=",sample_size,",init=",model_name,"_init,refresh=10,control=list(max_treedepth=",max_depth,",adapt_delta=",adapt_delta,"))"),
                 paste0("D_",model_name,"$inference = 1"),
                 paste0("S_",model_name," = stan('",model_name,".stan',data=D_",model_name,",chains=",chains,",iter=",sample_size,",init=",model_name,"_init,refresh=10,control=list(max_treedepth=",max_depth,",adapt_delta=",adapt_delta,"))"),
                 paste0("save(SIM_",model_name,",S_",model_name,",D_",model_name,",file='samples_",model_name,"_",gsub(" |:","_",Sys.time()),".Rdata')")
                 )
    writeLines(tt_Rfile,paste0("run_",model_name,".R"))
  }
}




bashfile_rdump_1C = function(model_name,data_file,warmup=500,iter=100,adapt_delta=0.8,max_depth=12,init=0,type="cmdstan",timelimit=96,chains=10) {
  if(type=="cmdstan") {
    data_file$inference=0
    with(data_file,stan_rdump(ls(data_file),file=paste0("onecountry/data_SIM_",model_name,".R")))
    data_file$inference=1
    with(data_file,stan_rdump(ls(data_file),file=paste0("onecountry/data_S_",model_name,".R")))
    if(length(init)!=1) {
      with(init,stan_rdump(ls(init),file=paste0("onecountry/init_",model_name,".R")))
      init = paste0("init_",model_name,".R")
    }
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=",chains),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_SIM_",model_name,".R output file=SIM_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10"),
                  paste0("./",model_name," sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,".R output file=S_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
    )
    writeLines(tt_launch,paste0("onecountry/sb_",model_name,".sh"))
  }
  if(type=="rstan") {
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  "module load R",
                  paste0("Rscript run_",model_name,".R"))
    writeLines(tt_launch,paste0("sbr_",model_name,".sh"))
    tt_Rfile = c("library(rstan);options(mc.cores = parallel::detectCores())",
                 paste0("load('data_",model_name,".Rdata')"),
                 paste0("D_",model_name,"$inference = 0"),
                 paste0("SIM_",model_name," = stan('",model_name,".stan',data=D_",model_name,",chains=",chains,",iter=",sample_size,",init=",model_name,"_init,refresh=10,control=list(max_treedepth=",max_depth,",adapt_delta=",adapt_delta,"))"),
                 paste0("D_",model_name,"$inference = 1"),
                 paste0("S_",model_name," = stan('",model_name,".stan',data=D_",model_name,",chains=",chains,",iter=",sample_size,",init=",model_name,"_init,refresh=10,control=list(max_treedepth=",max_depth,",adapt_delta=",adapt_delta,"))"),
                 paste0("save(SIM_",model_name,",S_",model_name,",D_",model_name,",file='samples_",model_name,"_",gsub(" |:","_",Sys.time()),".Rdata')")
    )
    writeLines(tt_Rfile,paste0("run_",model_name,".R"))
  }
}


bashfile_rdump_M7u = function(model_name,data_file,warmup=500,iter=100,adapt_delta=0.8,max_depth=12,init=0,type="cmdstan",timelimit=96,chains=10) {
    data_file$inference=1
    with(data_file,stan_rdump(ls(data_file),file=paste0("models/data_S_",model_name,".R")))
    if(length(init)!=1) {
      with(init,stan_rdump(ls(init),file=paste0("models/init_",model_name,".R")))
      init = paste0("init_",model_name,".R")
    }
    tt_launch = c("#!/bin/bash",
                  paste0("#SBATCH --job-name='",model_name,"'"),
                  "#SBATCH --partition=all",
                  paste0("#SBATCH --array=1-",chains),
                  paste0("#SBATCH --cpus-per-task=",data_file$G),
                  paste0("#SBATCH --time=",timelimit,":00:00"),
                  "",
                  "module load Boost/1.66.0-foss-2018a",
                  paste0("./M7 sample num_warmup=",warmup," num_samples=",iter," \\"),
                  paste0("      adapt delta=",adapt_delta," algorithm=hmc engine=nuts max_depth=",max_depth," init=",init," \\"),
                  paste0("      data file=data_S_",model_name,".R output file=S_",model_name,"_",gsub(" |:","-",Sys.time()),"_${SLURM_ARRAY_TASK_ID}.csv refresh=10")
    )
    writeLines(tt_launch,paste0("models/sb_",model_name,".sh"))
}

# plots
plot_indicators = function(samples=NULL,
                           data.list,
                           key) {
  
  G = data.list$G
  S = data.list$S
  L = data.list$L
  start = data.list$starting_year
  comp=c("A_pred","B_pred","C_pred","D_pred")
  comp2 = c("Prevalence","On ART","Mortality","Population")
  y = rstan::extract(samples,pars="y")[[1]]
  A_pred = rstan::extract(samples,pars="A_pred")[[1]]
  B_pred = rstan::extract(samples,pars="B_pred")[[1]]
  C_pred = rstan::extract(samples,pars="C_pred")[[1]]
  D_pred = rstan::extract(samples,pars="D_pred")[[1]]
  
  prop = list(A_prop = A_pred / D_pred,
              B_prop = B_pred / A_pred,
              C_prop = C_pred / A_pred,
              D_prop = D_pred / 1e6)
  
  pred = NULL
  for(i in 1:4) {
    for(j in 1:G) {
      pred = rbind(pred,as.data.frame(t(apply(prop[[i]][,j,],2,qsum))) %>%
                     tbl_df() %>%
                     mutate(CountryID=j,
                            comp=comp[i],
                            comp2=comp2[i],
                            year=1:L+start))
    }
  }
  pred = pred %>%
    left_join(key) %>%
    bind_cols(obs=c(as.vector(t(data.list$A_est/data.list$D_est)),
                    as.vector(t(data.list$B_est/data.list$A_est)),
                    as.vector(t(data.list$C_est/data.list$A_est)),
                    as.vector(t(data.list$D_est/1e6)))) 
  
  g_indic =
    ggplot(pred) +
    annotate("point",x=2000,y=0,colour="transparent") +
    geom_point(aes(x=year,y=obs),colour="grey30",size=.4,shape=21,fill="white") +
    geom_ribbon(aes(x=year,ymin=`2.5%`,ymax=`97.5%`,fill=comp2),alpha=.2) +
    geom_line(aes(x=year,y=`50%`,colour=comp2),size=.6) +
    facet_grid(comp~Country.Name2,
               scales = "free_y",
               labeller=labeller(comp=c("A_pred"="Prevalence (prop.)",
                                        "B_pred"="ART (prop.)",
                                        "C_pred"="Mortality (prop.)",
                                        "D_pred"="Population (mil.)"))) +
    scale_x_continuous(breaks=c(2000,2008,2016),labels=c("'00","'08","'16"),limits=c(1999,2018)) +
    scale_fill_brewer(palette="Spectral",guide=FALSE) +
    scale_colour_brewer(palette="Spectral",guide=FALSE) +
    labs(x="Year",y=NULL) 
  return(g_indic)
}


plot_indicators2 = function(samples=NULL,
                           data.list,
                           key,
                           grey=FALSE) {
  
  G = data.list$G
  S = data.list$S
  L = data.list$L
  start = data.list$starting_year
  comp=c("A_pred","B_pred","C_pred","D_pred")
  comp2 = c("Prevalence","On ART","Mortality","Population")
  A_pred = rstan::extract(samples,pars="A_pred")[[1]]
  B_pred = rstan::extract(samples,pars="B_pred")[[1]]
  C_pred = rstan::extract(samples,pars="C_pred")[[1]]
  D_pred = rstan::extract(samples,pars="D_pred")[[1]]
  
  prop = list(A_prop = A_pred ,
              B_prop = B_pred ,
              C_prop = C_pred ,
              D_prop = D_pred )
  
  pred = NULL
  for(i in 1:4) {
    for(j in 1:G) {
      pred = rbind(pred,as.data.frame(t(apply(prop[[i]][,j,],2,qsum))) %>%
                     tbl_df() %>%
                     mutate(CountryID=j,
                            comp=comp[i],
                            comp2=comp2[i],
                            year=1:L+start))
    }
  }
  pred = pred %>%
    left_join(key) %>%
    bind_cols(obs=c(as.vector(t(data.list$A_est)),
                    as.vector(t(data.list$B_est)),
                    as.vector(t(data.list$C_est)),
                    as.vector(t(data.list$D_est)))) 
  
  g_indic =
    ggplot(pred) +
    annotate(geom = "point",x=2000,y=0,colour="white",alpha=0) +
    geom_point(aes(x=year,y=obs),colour="grey30",size=.4,shape=21,fill="white") +
    geom_ribbon(aes(x=year,ymin=`2.5%`,ymax=`97.5%`,fill=comp2),alpha=.2) +
    geom_line(aes(x=year,y=`50%`,colour=comp2),size=.6) +
    facet_wrap(~comp,
               scales = "free_y",
               nrow=1,
               labeller=labeller(comp=c("A_pred"="Prevalence (Mil.)",
                                        "B_pred"="ART (Mil.)",
                                        "C_pred"="Mortality (Mil.)",
                                        "D_pred"="Population (Mil.)"))) +
    scale_x_continuous(breaks=c(2000,2008,2016),labels=c("'00","'08","'16"),limits=c(1999,2018)) +
    scale_y_continuous(labels=function(x) paste0(x/1e6,"M")) +
    labs(x="Year",y=NULL)  
  if(!grey) {
    g_indic = g_indic   +  
      scale_fill_brewer(palette="Spectral",guide=FALSE) +
      scale_colour_brewer(palette="Spectral",guide=FALSE) 
  } else {
    g_indic = g_indic   +  
      scale_fill_manual(values=rep("grey10",4),guide=FALSE) +
      scale_colour_manual(values=rep("grey70",4),guide=FALSE)
  }
  return(g_indic)
}



plot_pdr = function(sim.samples=NULL,
                           samples=NULL,
                           data.list,
                           showprior=TRUE,
                           ylim=NULL,
                           key,
                           ...) {
  G = data.list$G
  S = data.list$S
  L = data.list$L
  start = data.list$starting_year
  comp=c("F_output")
  comp2 = c("Pretreatment drug resistance")
  if(dim(key)[[1]]==1) key$CountryID = 1
  
  if(!is.null(sim.samples)) {
    pred = summary(sim.samples,pars=comp)[[1]] %>%
      as.data.frame(.) %>%
      rownames_to_column() %>%
      tbl_df() %>%
      mutate(year=start+rep(1:L,G),
             CountryID=rep(1:G,each=L),
             comp=rep(comp,each=G*L),
             comp2=rep(comp2,each=G*L)) %>%
      left_join(key)
    cc = "gray30"
  }
  
  if(!is.null(samples)) {
    pred = summary(samples,pars=comp)[[1]] %>%
      as.data.frame(.) %>%
      rownames_to_column() %>%
      tbl_df() %>%
      mutate(year=start+rep(1:L,G),
             CountryID=rep(1:G,each=L),
             comp=rep(comp,each=G*L),
             comp2=rep(comp2,each=G*L)) %>%
      left_join(key)
    cc = "purple"
  }
  
  obs = data.frame(k=data.list$F_k,
                   n=data.list$F_n,
                   CountryID=data.list$F_country,
                   year=data.list$F_t+start) %>%
    tbl_df() %>% 
    mutate(p=k/n,
           pinf=qbeta(0.025, k, n-k+1),
           psup=qbeta(0.975, k+1, n-k)) %>%
    left_join(key) %>%
    group_by(CountryID,year) %>%
      mutate(rank=row_number(),
             year2=year+(rank-1)*.2) %>%
    ungroup()
  g =
    ggplot() +
    geom_ribbon(data=pred,aes(x=year,ymin=`2.5%`,ymax=`97.5%`),fill=cc,alpha=.2) +
    geom_pointrange(data=obs,aes(x=year2,y=p,ymin=pinf,ymax=psup),colour="grey30",size=.2,shape=21,fill="white",stroke=.7) +
    facet_wrap(~Country.Name2,
               ncol=3) +
    scale_x_continuous(breaks=c(2002,2006,2010,2014,2018),labels=c("'02","'06","'10","'14","'18")) +
    scale_y_continuous(labels=scales::percent) +
    labs(x="Year",y=comp2) 
  if(!is.null(ylim)) g = g + coord_cartesian(ylim=c(0,ylim))
  if(!is.null(samples)) g = g + geom_line(data=pred,aes(x=year,y=`50%`),alpha=.6,size=.4)
  return(g)
}



plot_pdr2 = function(sim.samples=NULL,
                    samples=NULL,
                    data.list,
                    showprior=TRUE,
                    ylim=NULL,
                    key,
                    ...) {
  G = data.list$G
  S = data.list$S
  L = data.list$L
  start = data.list$starting_year
  comp=c("F_output")
  comp2 = c("Pretreatment drug resistance")
  if(dim(key)[[1]]==1) key$CountryID = 1
  
  if(!is.null(sim.samples)) {
    pred = summary(sim.samples,pars=comp)[[1]] %>%
      as.data.frame(.) %>%
      rownames_to_column() %>%
      tbl_df() %>%
      mutate(year=start+rep(1:L,G),
             CountryID=rep(1:G,each=L),
             comp=rep(comp,each=G*L),
             comp2=rep(comp2,each=G*L)) %>%
      left_join(key)
    cc = "gray30"
  }
  
  if(!is.null(samples)) {
    pred = summary(samples,pars=comp)[[1]] %>%
      as.data.frame(.) %>%
      rownames_to_column() %>%
      tbl_df() %>%
      mutate(year=start+rep(1:L,G),
             CountryID=rep(1:G,each=L),
             comp=rep(comp,each=G*L),
             comp2=rep(comp2,each=G*L)) %>%
      left_join(key)
    cc = "seagreen"
  }
  
  obs = data.frame(k=data.list$F_k,
                   n=data.list$F_n,
                   CountryID=data.list$F_country,
                   year=data.list$F_t+start) %>%
    tbl_df() %>% 
    mutate(p=k/n,
           pinf=qbeta(0.025, k, n-k+1),
           psup=qbeta(0.975, k+1, n-k)) %>%
    left_join(key) %>%
    group_by(CountryID,year) %>%
    mutate(rank=row_number(),
           year2=year+(rank-1)*.2) %>%
    ungroup()
  g =
    ggplot() +
    annotate("point",x=start,y=0,col="white") +
    geom_ribbon(data=pred,aes(x=year,ymin=`2.5%`,ymax=`97.5%`),fill=cc,alpha=.2) +
    geom_pointrange(data=obs,aes(x=year2,y=p,ymin=pinf,ymax=psup),colour="grey30",size=.2,shape=21,fill="white",stroke=.7) +
    scale_x_continuous(breaks=c(2000,2008,2016),labels=c("'00","'08","'16")) +
    scale_y_continuous(labels=scales::percent) +
    labs(x="Year",y=comp2) 
  if(!is.null(ylim)) g = g + coord_cartesian(ylim=c(0,ylim))
  if(!is.null(samples)) g = g + geom_line(data=pred,aes(x=year,y=`50%`),colour=cc,size=.4)
  return(g)
}

plot_compartments = function(samples=NULL,
                           data.list,
                           key) {
  
  G = data.list$G
  S = data.list$S
  L = data.list$L
  C = data.list$C
  start = data.list$starting_year
  y = rstan::extract(samples,pars="y")[[1]]

  pred = NULL
  for(i in 1:C) {
    for(j in 1:G) {
      pred = rbind(pred,as.data.frame(t(apply(y[,j,,i],2,qsum))) %>%
                     tbl_df() %>%
                     mutate(CountryID=j,
                            comp=i,
                            year=1:S+2000))
    }
  }
  pred = left_join(pred,key)

  g_indic =
    ggplot(pred) +
    geom_ribbon(aes(x=year,ymin=`2.5%`,ymax=`97.5%`,fill=as.factor(comp)),alpha=.3) +
    geom_line(aes(x=year,y=`50%`)) +
    facet_grid(comp~Country.Name2,scales="free") +
    scale_x_continuous(breaks=c(2000,2008,2016),labels=c("'00","'08","'16"),limits=c(2000,2016)) +
    scale_fill_brewer(palette="Spectral",guide=FALSE) +
    scale_colour_brewer(palette="Spectral",guide=FALSE) +
    labs(x="Year",y=NULL)
  return(g_indic)
}

check_compartments = function(samples=NULL,
                             data.list,
                             key) {
  
  G = data.list$G
  S = data.list$S
  L = data.list$L
  C = data.list$C
  start = data.list$starting_year
  y = rstan::extract(samples,pars="y")[[1]][,,1:16,]
  lapply(1:10,function(x,...) sum(is.na(y[,,1,x])))
  (miss_y = which(is.na(y[,,1,1])))
  output = rstan::extract(samples,pars=c("A_output","B_output","C_output","D_output","E_output"))
  lapply(output,function(x) sum(is.na(x[,,1])))
  (miss_output = which(is.na(output[[1]][,,1])))
  cat("same NAs:",all.equal(miss_y,miss_output))
  pars = rstan::extract(samples,pars=c("beta","tau","nu","xi","eta","delta","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota"))
  lapply(pars,function(x) sum(is.na(x)))
  (miss_pars = which(is.na(pars[[1]])))
  cat("same NAs:",all.equal(miss_y,miss_pars))
  rawpars = rstan::extract(samples,pars=c("beta_raw","tau_raw","nu_raw","xi_raw","eta_raw","delta_raw","omega_raw","iota_raw"))
  lapply(rawpars,function(x) sum(is.na(x)))
  
  for(i in 1:8) {
    plot(density(rawpars[[i]]),main=i)
    points(density(rawpars[[i]][miss_pars]),col="red",type="l")
  }
  
  checkpars = list()
  checkpars[[1]] = rawpars[[1]] / data.list$p_beta
  checkpars[[2]] = rawpars[[2]] / data.list$p_tau
  checkpars[[3]] = rawpars[[3]] + 0.5
  checkpars[[4]] = data.list$t1 + rawpars[[4]] / data.list$p_nu
  checkpars[[5]] = rawpars[[5]] / data.list$p_eta
  checkpars[[6]] = rawpars[[6]] / data.list$p_delta
  checkpars[[7]] = lognormal_repar(rawpars[[7]],c(data.list$p_mu_omega,data.list$p_sigma_omega))
  checkpars[[8]] = inv.logit(rawpars[[8]])
  
  for(i in 1:8) {
    plot(density(checkpars[[i]]),main=i)
    points(density(checkpars[[i]][miss_pars]),col="red",type="l")
  }
}
