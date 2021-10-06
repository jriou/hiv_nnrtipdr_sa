## Set-up -----------------------------------------------------------------------

load("all_data.Rdata")
source("model-dev_setup.R")

colsurvey = "darkorange"
colindic = "maroon"

starting_year = 1999
pdr_data = filter(updated_data,year>=starting_year)
country_data = filter(country_data,year>=starting_year)


### Model M1 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M1 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=5,D=2,E=2,G=9)
D_M1$p_beta = 5
D_M1$p_tau = 5
D_M1$p_sigma = 1

SIM_M1 = stan("models/M1.stan",data=D_M1,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M1",data_file=D_M1,warmup=1000,iter=1000,adapt_delta=0.9,max_depth=15,init=1,type="cmdstan",timelimit=24,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M1.stan \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M1.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M1.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M1.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M1
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M1.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M1 = read_stan_csv(dir(".",pattern = 'SIM_M1_2019-07-09-10-52-41_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M1)
S_M1 = read_stan_csv(dir(".",pattern = 'S_M1_2019-07-09-10-52-41_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M1)
D_S_M1 = read_rdump("data_S_M1.R")
D_SIM_M1 = read_rdump("data_SIM_M1.R")
save(SIM_M1,D_SIM_M1,S_M1,D_S_M1,file="M1_2019-07-09.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M1_2019-07-09.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M1_2019-07-09.Rdata")


check_hmc_diagnostics(SIM_M1)
ppM1 = c("beta","tau","sigma")
print(SIM_M1,pars=ppM1,digits_summary = 4)
plot_indicators(samples=SIM_M1,data.list=D_SIM_M1,key=country_key)


check_hmc_diagnostics(S_M1)
ppM1 = c("beta","tau")
print(S_M1,pars=ppM1,digits_summary = 4)
plot_indicators(samples=S_M1,data.list=D_S_M1,key=country_key)
loo(S_M1)
stan_trace(S_M1,c("beta","tau"))





### Model M2 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M2 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=5,D=4,E=2,G=9)
D_M2$p_beta = 5
D_M2$p_tau = 5
D_M2$p_xi = 2
D_M2$p_nu = .5
D_M2$p_sigma = 1

SIM_M2 = stan("models/M2.stan",data=D_M2,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M2",data_file=D_M2,warmup=1000,iter=1000,adapt_delta=0.9,max_depth=15,init=1,type="cmdstan",timelimit=24,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M2.stan \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M2.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M2.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M2.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M2
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M2.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M2 = read_stan_csv(dir(".",pattern = 'SIM_M2_2019-11-12-12-21-41_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M2)
S_M2 = read_stan_csv(dir(".",pattern = 'S_M2_2019-11-12-12-21-41_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M2)
D_S_M2 = read_rdump("data_S_M2.R")
D_SIM_M2 = read_rdump("data_SIM_M2.R")
save(SIM_M2,D_SIM_M2,S_M2,D_S_M2,file="M2_2019-11-12.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M2_2019-11-12.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M2_2019-11-12.Rdata")


check_hmc_diagnostics(SIM_M2)
ppM2 = c("beta","tau")
print(SIM_M2,pars=ppM2,digits_summary = 4)
plot_indicators(samples=SIM_M2,data.list=D_SIM_M2,key=country_key)

print(S_M2,pars="A_output")

check_hmc_diagnostics(S_M2)
ppM2 = c("beta","tau")
print(S_M2,pars=ppM2,digits_summary = 4)
plot_indicators(samples=S_M2,data.list=D_S_M2,key=country_key)
plot_compartments(samples=S_M2,data.list=D_S_M2,key=country_key)
loo(S_M2)




### Model M3 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M3 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=5,D=6,E=4,G=9)
D_M3$p_beta = 5
D_M3$p_tau = 5
D_M3$p_xi = 2
D_M3$p_nu = .5
D_M3$p_eta = 10
D_M3$p_delta = 10
D_M3$p_sigma = 1


SIM_M3 = stan("models/M3.stan",data=D_M3,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M3",data_file=D_M3,warmup=1000,iter=1000,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=24,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M3.stan \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M3.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M3.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M3.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M3
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M3.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
  SIM_M3 = read_stan_csv(dir(".",pattern = 'SIM_M3_2019-11-12-16-24-13_[[:digit:]]+.csv'))
  check_hmc_diagnostics(SIM_M3)
  S_M3 = read_stan_csv(dir(".",pattern = 'S_M3_2019-11-12-16-24-13_[[:digit:]]+.csv'))
  check_hmc_diagnostics(S_M3)
  print(S_M3,pars="delta")
  D_S_M3 = read_rdump("data_S_M3.R")
  D_SIM_M3 = read_rdump("data_SIM_M3.R")
  save(SIM_M3,D_SIM_M3,S_M3,D_S_M3,file="M3_2019-11-12.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M3_2019-11-12.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M3_2019-11-12.Rdata")


check_hmc_diagnostics(SIM_M3)
ppM3 = c("beta","tau","nu","xi","delta","eta")
print(SIM_M3,pars=ppM3,digits_summary = 4)
plot_indicators(samples=SIM_M3,data.list=D_S_M3,key=country_key)


check_hmc_diagnostics(S_M3)
ppM3 = c("beta","tau","nu","xi","delta","eta")
print(S_M3,pars=ppM3,digits_summary = 4)
plot_indicators(samples=S_M3,data.list=D_S_M3,key=country_key)
stan_plot(S_M3,pars="beta")
stan_plot(S_M3,pars="tau")
stan_plot(S_M3,pars="nu")
stan_plot(S_M3,pars="xi")
stan_plot(S_M3,pars="eta")
stan_plot(S_M3,pars="delta")
stan_plot(S_M3,pars="sigma")
stan_trace(S_M3,pars="eta")
stan_dens(S_M3,pars="eta",separate_chains =TRUE)



### Model M4 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M4 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=10,D=8,E=4,G=9)
D_M4$p_beta = 5
D_M4$p_tau = 5
D_M4$p_xi = 2
D_M4$p_nu = .5
D_M4$p_eta = 5
D_M4$p_delta = 5
D_M4$p_sigma = 1
D_M4$p_mu_omega = 5
D_M4$p_sigma_omega = 20
D_M4$p_mu_iota = c(1,9)
D_M4$p_sigma_iota = 5

SIM_M4 = stan("models/M4.stan",data=D_M4,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M4",data_file=D_M4,warmup=1000,iter=1000,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=24,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M4.stan \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M4.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M4.R \
    /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M4.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M4
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M4.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M4 = read_stan_csv(dir(".",pattern = 'SIM_M4_2019-11-13-16-36-33_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M4)
S_M4 = read_stan_csv(dir(".",pattern = 'S_M4_2019-11-13-16-36-33_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M4)
D_S_M4 = read_rdump("data_S_M4.R")
D_SIM_M4 = read_rdump("data_SIM_M4.R")
save(SIM_M4,D_SIM_M4,S_M4,D_S_M4,file="M4_2019-11-13.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M4_2019-11-13.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M4_2019-11-13.Rdata")


check_hmc_diagnostics(SIM_M4)
ppM4 = c("beta","tau","nu","xi","eta","delta","omega","iota")
print(SIM_M4,pars=ppM4,digits_summary = 4)
plot_indicators(samples=SIM_M4,data.list=D_M4,key=country_key)
plot_pdr(sim.samples=SIM_M4,data.list=D_M4,key=country_key)
plot_compartments(SIM_M4,data.list=D_SIM_M4,key=country_key)

check_hmc_diagnostics(S_M4)
ppM4 = c("beta","tau","nu","xi","eta","delta","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M4,pars=ppM4,digits_summary = 4)

print(S_M4,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 2)
print(S_M4,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)

plot_indicators(samples=S_M4,data.list=D_S_M4,key=country_key)
ggsave(file="figures/post_indicators_M4.pdf",width=8,height=5)
plot_pdr(samples=S_M4,data.list=D_S_M4,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M4.pdf",width=5.5,height=4)
loo(S_M4)
loo(S_M4,pars="log_lik2")
plot_compartments(samples=S_M4,data.list=D_S_M4,key=country_key)



summary(S_M4,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  # filter(CountryID != c(3,8)) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`75%`,ymin=`25%`))

summary(S_M4,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))




### Model M5 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M5 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=10,D=8,E=4,G=9)
D_M5$p_beta = 5
D_M5$p_tau = 5
D_M5$p_xi = 2
D_M5$p_nu = .5
D_M5$p_eta = 5
D_M5$p_delta = 5
D_M5$p_sigma = 1
D_M5$p_mu_omega = 5
D_M5$p_sigma_omega = 20
D_M5$p_mu_iota = c(1,9)
D_M5$p_sigma_iota = 5
D_M5$P = 8
D_M5$covariables = country_char_matrix
D_M5$p_pi = 1

SIM_M5 = stan("models/M5.stan",data=D_M5,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M5",data_file=D_M5,warmup=1000,iter=1000,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=48,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M5.stan \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M5.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M5.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M5.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M5
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M5.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M5 = read_stan_csv(dir(".",pattern = 'SIM_M5_2019-11-12-17-18-44_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M5)
S_M5 = read_stan_csv(dir(".",pattern = 'S_M5_2019-11-13-17-18-44_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M5)
D_S_M5 = read_rdump("data_S_M5.R")
D_SIM_M5 = read_rdump("data_SIM_M5.R")
save(S_M5,D_S_M5,file="M5_2019-11-13.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M5_2019-11-13.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M5_2019-11-13.Rdata")


check_hmc_diagnostics(SIM_M5)
ppM5 = c("beta","tau","nu","xi","eta","delta","omega","iota","pi")
print(SIM_M5,pars=ppM5,digits_summary = 4)
plot_indicators(samples=SIM_M5,data.list=D_M5,key=country_key)
plot_pdr(sim.samples=SIM_M5,data.list=D_M5,key=country_key)
plot_compartments(SIM_M5,data.list=D_SIM_M5,key=country_key)

y = extract(SIM_M5,pars="y")[[1]]
sum(y<0,na.rm=T)
y[which(y<0)]
sum(y[,,12:19,2]<0,na.rm=T)
pp = extract(SIM_M5,pars=ppM5)

check_hmc_diagnostics(S_M5)
ppM5 = c("beta","tau","nu","xi","eta","delta","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M5,pars=ppM5,digits_summary = 4)

print(S_M5,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 2)
print(S_M5,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)
print(S_M5,pars=c("pi"),digits_summary = 3)

plot_indicators(samples=S_M5,data.list=D_S_M5,key=country_key)
ggsave(file="figures/post_indicators_M5.pdf",width=8,height=5)
plot_pdr(samples=S_M5,data.list=D_S_M5,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M5.pdf",width=5.5,height=4)
loo(S_M5)
plot_compartments(samples=S_M5,data.list=D_S_M5,key=country_key)



summary(S_M5,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  # filter(CountryID != c(3,8)) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`75%`,ymin=`25%`))

summary(S_M5,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))




### Model M6 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M6 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=11,D=9,E=4,G=9)
D_M6$p_beta = 5
D_M6$p_tau = 5
D_M6$p_xi = 2
D_M6$p_nu = .5
D_M6$p_eta = 5
D_M6$p_delta = 5
D_M6$p_kappa = 5
D_M6$p_sigma = 1
D_M6$p_mu_omega = 5
D_M6$p_sigma_omega = 20
D_M6$p_mu_iota = c(1,9)
D_M6$p_sigma_iota = 5

SIM_M6 = stan("models/M6.stan",data=D_M6,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M6",data_file=D_M6,warmup=500,iter=500,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=96,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M6.stan \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M6.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M6.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M6.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M6
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M6.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
# SIM_M6 = read_stan_csv(dir(".",pattern = 'SIM_M6_2020-06-10-13-49-51_[[:digit:]]+.csv'))
# check_hmc_diagnostics(SIM_M6)
S_M6 = read_stan_csv(dir(".",pattern = 'S_M6_2020-06-10-13-49-51_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M6)
D_S_M6 = read_rdump("data_S_M6.R")
# D_SIM_M6 = read_rdump("data_SIM_M6.R")
save(S_M6,D_S_M6,file="M6_2020-06-10.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
system("scp UBELIX:projects/nnrti_sa3/models/M6_2020-06-10.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.")
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M6_2020-06-10.Rdata")


check_hmc_diagnostics(SIM_M6)
ppM6 = c("beta","tau","nu","xi","eta","delta","omega","iota")
print(SIM_M6,pars=ppM6,digits_summary = 4)
plot_indicators(samples=SIM_M6,data.list=D_M6,key=country_key)
plot_pdr(sim.samples=SIM_M6,data.list=D_M6,key=country_key)
plot_compartments(SIM_M6,data.list=D_SIM_M6,key=country_key)

check_hmc_diagnostics(S_M6)
ppM6 = c("beta","tau","nu","xi","eta","delta","kappa","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M6,pars=ppM6,digits_summary = 4)

print(S_M6,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 4)
print(S_M6,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)

plot_indicators(samples=S_M6,data.list=D_S_M6,key=country_key)
ggsave(file="figures/post_indicators_M6.pdf",width=8,height=5)
plot_pdr(samples=S_M6,data.list=D_S_M6,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M6.pdf",width=5.5,height=4)
loo(S_M6)
loo(S_M6,pars="log_lik2")
plot_compartments(samples=S_M6,data.list=D_S_M6,key=country_key)


compare(loo(S_M4,pars="log_lik"),loo(S_M6,pars="log_lik"))
compare(loo(S_M4,pars="log_lik2"),loo(S_M6,pars="log_lik2"))


summary(S_M6,pars="kappa")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  # filter(CountryID != c(3,8)) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`75%`,ymin=`25%`))

summary(S_M6,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))


summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  arrange(mean)

summary(S_M6,pars="nu")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  arrange(`50%`)
  

### Model M7uu ---------------------------------------------------------------------------------------------------
# transform data in list form
for(i in 1:9) {
    
  D_M7u = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                         S=20,C=11,D=9,E=4,G=9)
  D_M7u$p_beta = 5
  D_M7u$p_tau = 5
  D_M7u$p_xi = 2
  D_M7u$p_nu = .5
  D_M7u$p_eta = 5
  D_M7u$p_delta = 5
  D_M7u$p_kappa = 5
  D_M7u$p_sigma = 1
  D_M7u$p_mu_omega = 5
  D_M7u$p_sigma_omega = 20
  D_M7u$p_mu_iota = c(1,9)
  D_M7u$p_sigma_iota = 5
  D_M7u$P = 1
  D_M7u$covariables = structure(country_char_matrix[,i],.Dim=c(9,1))
  D_M7u$p_pi = 1
  D_M7u$cov_name = i
  
  # inference
  # >>>>>>>>>>>>>>>>> prepare for cluster:
  bashfile_rdump_M7u(model_name=paste0("M7u",i),data_file=D_M7u,warmup=500,iter=500,adapt_delta=0.8,max_depth=20,init=0,type="cmdstan",timelimit=24,chains=4)
}


# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M7u*.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M7u*.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M7
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M7u.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
for(i in 6:9) {
  assign(paste0("S_M7u",i),read_stan_csv(dir(".",pattern = paste0('S_M7u',i,'_2019-11-29-18-02-36_[[:digit:]]+.csv'))))
  assign(paste0("D_S_M7u",i),read_rdump(paste0("data_S_M7u",i,".R")))
  print(get(paste0("S_M7u",i)),pars="pi")
  cat(i)
}
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M7u*_2019-11-28.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M7u_2019-11-28.Rdata")


check_hmc_diagnostics(S_M7u)
ppM7u = c("beta","tau","nu","xi","eta","delta","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M7u,pars=ppM7u,digits_summary = 4)

print(S_M7u,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 2)
print(S_M7u,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)
print(S_M7u,pars=c("pi"),digits_summary = 3)

sig = extract(S_M7u,pars="sigma_omega")[[1]] / extract(S_M6,pars="sigma_omega")[[1]]
qsum(sig)

plot_indicators(samples=S_M7u,data.list=D_S_M7u,key=country_key)
ggsave(file="figures/post_indicators_M7u.pdf",width=8,height=5)
plot_pdr(samples=S_M7u,data.list=D_S_M7u,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M7u.pdf",width=5.5,height=4)
loo(S_M7u)
plot_compartments(samples=S_M7u,data.list=D_S_M7u,key=country_key)

### Model M7 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M7 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=11,D=9,E=4,G=9)
D_M7$p_beta = 5
D_M7$p_tau = 5
D_M7$p_xi = 2
D_M7$p_nu = .5
D_M7$p_eta = 5
D_M7$p_delta = 5
D_M7$p_kappa = 5
D_M7$p_sigma = 1
D_M7$p_mu_omega = 5
D_M7$p_sigma_omega = 20
D_M7$p_mu_iota = c(1,9)
D_M7$p_sigma_iota = 5
D_M7$P = 9
D_M7$covariables = country_char_matrix
D_M7$p_pi = 1

SIM_M7 = stan("models/M7.stan",data=D_M7,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M7",data_file=D_M7,warmup=1000,iter=1000,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=48,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M7.stan \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M7.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M7.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M7.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M7
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M7.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M7 = read_stan_csv(dir(".",pattern = 'SIM_M7_2019-11-28-15-05-13_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M7)
S_M7 = read_stan_csv(dir(".",pattern = 'S_M7_2019-11-28-15-05-13_[[:digit:]]+.csv'))
check_hmc_diagnostics(S_M7)
D_S_M7 = read_rdump("data_S_M7.R")
D_SIM_M7 = read_rdump("data_SIM_M7.R")
save(S_M7,D_S_M7,file="M7_2019-11-28.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M7_2019-11-28.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M7_2019-11-28.Rdata")

  
check_hmc_diagnostics(SIM_M7)
ppM7 = c("beta","tau","nu","xi","eta","delta","omega","iota","pi")
print(SIM_M7,pars=ppM7,digits_summary = 4)
plot_indicators(samples=SIM_M7,data.list=D_M7,key=country_key)
plot_pdr(sim.samples=SIM_M7,data.list=D_M7,key=country_key)
plot_compartments(SIM_M7,data.list=D_SIM_M7,key=country_key)


check_hmc_diagnostics(S_M7)
ppM7 = c("beta","tau","nu","xi","eta","delta","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M7,pars=ppM7,digits_summary = 4)

print(S_M7,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 2)
print(S_M7,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)
print(S_M7,pars=c("pi"),digits_summary = 3)

sig = extract(S_M7,pars="sigma_omega")[[1]] / extract(S_M6,pars="sigma_omega")[[1]]
qsum(sig)

plot_indicators(samples=S_M7,data.list=D_S_M7,key=country_key)
ggsave(file="figures/post_indicators_M7.pdf",width=8,height=5)
plot_pdr(samples=S_M7,data.list=D_S_M7,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M7.pdf",width=5.5,height=4)
loo(S_M7)
plot_compartments(samples=S_M7,data.list=D_S_M7,key=country_key)

compare(loo(S_M6,"log_lik"),loo(S_M7,"log_lik"))
compare(loo(S_M6,"log_lik2"),loo(S_M7,"log_lik2"))

summary(S_M7,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  # filter(CountryID != c(3,8)) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`75%`,ymin=`25%`))

summary(S_M7,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))





### Model M8 ---------------------------------------------------------------------------------------------------

# transform data in list form
D_M8 = generate_base_D(indic_data=country_data,survey_data=pdr_data,key=country_key,starting_year=starting_year,
                       S=20,C=11,D=8,E=4,G=9)
D_M8$p_beta = 5
D_M8$p_tau = 5
D_M8$p_xi = 2
D_M8$p_nu = .5
D_M8$p_eta = 5
D_M8$p_delta = 5
D_M8$p_kappa = 5
D_M8$p_sigma = 1
D_M8$p_mu_omega = 5
D_M8$p_sigma_omega = 20
D_M8$p_iota = c(1,9)

SIM_M8 = stan("models/M8.stan",data=D_M8,chains=1,iter=10)

# inference
# >>>>>>>>>>>>>>>>> prepare for cluster:
bashfile_rdump(model_name="M8",data_file=D_M8,warmup=1000,iter=1000,adapt_delta=0.95,max_depth=20,init=1,type="cmdstan",timelimit=24,chains=4)
# >>>>>>>>>>>>>>>>> export to cluster:
scp /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/M8.stan \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_SIM_M8.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/data_S_M8.R \
/home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/sb_M8.sh UBELIX:projects/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> compile:
make ../projects/nnrti_sa3/models/M8
# >>>>>>>>>>>>>>>>> launch:
sbatch sb_M8.sh
# >>>>>>>>>>>>>>>>> check & format data:
ml R
R
library(rstan)
SIM_M8 = read_stan_csv(dir(".",pattern = 'SIM_M8_2019-11-18-11-47-14_[[:digit:]]+.csv'))
check_hmc_diagnostics(SIM_M8)
S_M8 = read_stan_csv(dir(".",pattern = 'S_M8_2019-11-28-10-52-23_[[:digit:]]+.csv')[c(1,3:4)])
check_hmc_diagnostics(S_M8)
D_S_M8 = read_rdump("data_S_M8.R")
D_SIM_M8 = read_rdump("data_SIM_M8.R")
save(S_M8,D_S_M8,file="M8_2019-11-28.Rdata")
# >>>>>>>>>>>>>>>>> import data back:
scp UBELIX:projects/nnrti_sa3/models/M8_2019-11-28.Rdata /home/julien/Dropbox/Unibe/hiv_res/nnrti_sa3/models/.
# >>>>>>>>>>>>>>>>> load data:
l = load("models/M8_2019-11-28.Rdata")


check_hmc_diagnostics(SIM_M8)
ppM8 = c("beta","tau","nu","xi","eta","delta","omega","iota")
print(SIM_M8,pars=ppM8,digits_summary = 4)
plot_indicators(samples=SIM_M8,data.list=D_M8,key=country_key)
plot_pdr(sim.samples=SIM_M8,data.list=D_M8,key=country_key)
plot_compartments(SIM_M8,data.list=D_SIM_M8,key=country_key)

check_hmc_diagnostics(S_M8)
ppM8 = c("beta","tau","nu","xi","eta","delta","kappa","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M8,pars=ppM8,digits_summary = 4)

print(S_M8,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 4)
print(S_M8,pars=c("iota"),digits_summary = 3)

plot_indicators(samples=S_M8,data.list=D_S_M8,key=country_key)
ggsave(file="figures/post_indicators_M8.pdf",width=8,height=5)
plot_pdr(samples=S_M8,data.list=D_S_M8,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M8.pdf",width=5.5,height=4)
loo(S_M8)
loo(S_M8,pars="log_lik2")
plot_compartments(samples=S_M8,data.list=D_S_M8,key=country_key)


compare(loo(S_M6,pars="log_lik"),loo(S_M8,pars="log_lik"))
compare(loo(S_M6,pars="log_lik2"),loo(S_M8,pars="log_lik2"))


summary(S_M8,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  # filter(CountryID != c(3,8)) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`75%`,ymin=`25%`))

summary(S_M8,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))


summary(S_M8,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  arrange(mean)

summary(S_M8,pars="nu")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  arrange(`50%`)




