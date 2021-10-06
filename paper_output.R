## Set-up -----------------------------------------------------------------------

load("all_data.Rdata")
source("model-dev_setup.R")
library(xtable)

colsurvey = "darkorange"
colindic = "maroon"

starting_year = 1999
pdr_data = filter(updated_data,year>starting_year)
country_data = filter(country_data,year>=starting_year)
country_key$Country.Abb2 = c("BWA","LSO","MWI","MOZ","NAM","RSA","ESW","ZMB","ZWE")

## Model M6 -----------------------------------------------------------------------

l = load("models/M6_2020-06-10.Rdata")

check_hmc_diagnostics(S_M6)
ppM6 = c("beta","tau","nu","xi","eta","delta","kappa","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M6,pars=ppM6,digits_summary = 4)

print(S_M6,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 4)
print(S_M6,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)
print(S_M6,pars=c("kappa"),digits_summary = 3)

plot_pdr(samples=S_M6,data.list=D_S_M6,key=country_key,ylim = .3)


## Format relevant samples ----------------------------------------------

tt_tau = rstan::extract(S_M6,pars="tau")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","tau",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_beta = rstan::extract(S_M6,pars="beta")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","beta",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_xi = rstan::extract(S_M6,pars="xi")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","xi",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_nu = rstan::extract(S_M6,pars="nu")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","nu",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_omega = rstan::extract(S_M6,pars="omega")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","omega",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_iota = rstan::extract(S_M6,pars="iota")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","iota",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_kappa = rstan::extract(S_M6,pars="kappa")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","kappa",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_cum_art = rstan::extract(S_M6,pars="y")[[1]][,,18,9] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","cum_art",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
init_prev = country_key %>%
  mutate(init_prev = D_S_M6$init_prev)
tt_cum_inc = rstan::extract(S_M6,pars="y")[[1]][,,18,7] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","cum_inc",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  left_join(init_prev) %>%
  select(-pop) %>%
  left_join(country_key) %>%
  mutate(cum_inc=((cum_inc*pop)+init_prev)/pop) %>%
  select(CountryID,cum_inc) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_pred2018 = rstan::extract(S_M6,pars="F_output")[[1]][,,19] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","pred2018",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_delta = rstan::extract(S_M6,pars="delta")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","delta",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())
tt_cum_inc = rstan::extract(S_M6,pars="y")[[1]][,,18,7] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  gather("CountryID","cum_inc",1:9) %>%
  mutate(CountryID=as.numeric(gsub("V","",CountryID))) %>%
  left_join(init_prev) %>%
  left_join(country_key) %>%
  mutate(cum_inc=((cum_inc*pop)+init_prev)/pop) %>%
  select(CountryID,cum_inc) %>%
  group_by(CountryID) %>%
  mutate(it=row_number())

country_key$pop = D_S_M6$init_pop
tt = left_join(tt_tau,tt_nu) %>%
  left_join(tt_xi) %>%
  left_join(tt_beta) %>%
  left_join(tt_omega) %>%
  left_join(tt_iota) %>%
  left_join(tt_kappa) %>%
  left_join(tt_delta) %>%
  left_join(tt_cum_art) %>%
  left_join(tt_cum_inc) %>%
  left_join(tt_pred2018) %>%
  left_join(country_key) %>%
  mutate(argnt=omega/kappa)


## Text ----------------------------------------------------------------------------

omega = extract(S_M6,pars="omega")[[1]]
mu_omega = extract(S_M6,pars="mu_omega")[[1]]
kappa = extract(S_M6,pars="kappa")[[1]]
comp = omega
for(i in 1:9) {
  comp[,i] = omega[,i] / kappa[,i]
}

adjusted_growth = t(apply(comp,2,qsum)) %>%
  as.data.frame()  %>%
  mutate(CountryID=1:9,
         est=paste0(sprintf("%.2f",`50%`)," (",sprintf("%.2f",`2.5%`)," - ",sprintf("%.2f",`97.5%`),")")) %>%
  arrange(`50%`) %>%
  left_join(country_key) 
adjusted_growth
ggplot(adjusted_growth) +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`)) +
  scale_y_log10(expand=c(0,0)) +
  coord_cartesian(ylim=c(0.001,11))

summary(S_M6,pars="mu_omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  select(`50%`,`2.5%`,`97.5%`) %>%
  mutate(inv50=1-exp(-`50%`),inv25=1-exp(-`2.5%`),inv97=1-exp(-`97.5%`))

## omega
summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# iota
summary(S_M6,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# pred
summary(S_M6,pars="F_output")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=rep(1:9,each=D_S_M6$L),
         year=rep(1999+(1:D_S_M6$L),9)) %>%
  left_join(country_key) %>%
  filter(year==2018) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`) 

# iota in proportion of F_output

iota = extract(S_M6,pars="iota")[[1]]
p2018 = extract(S_M6,pars="F_output")[[1]][,,19]
comp = iota
for(i in 1:9) {
  comp[,i] = iota[,i] / p2018[,i]
}

prop_iota = t(apply(comp,2,qsum)) %>%
  as.data.frame()  %>%
  mutate(CountryID=1:9,
         est=paste0(sprintf("%.2f",`50%`)," (",sprintf("%.2f",`2.5%`)," - ",sprintf("%.2f",`97.5%`),")"),
         est2=paste0(sprintf("%.0f",`50%`*100),"% (",sprintf("%.0f",`2.5%`*100)," - ",sprintf("%.0f",`97.5%`*100),")")) %>%
  arrange(`50%`) %>%
  left_join(country_key) 
prop_iota
ggplot(prop_iota) +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymax=`97.5%`,ymin=`2.5%`))


# tau
summary(S_M6,pars="tau")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# nu
summary(S_M6,pars="nu")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# kappa
summary(S_M6,pars="kappa")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# beta
summary(S_M6,pars="beta")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)

# beta
summary(S_M6,pars="delta")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name2,`50%`,`2.5%`,`97.5%`) %>%
  arrange(`50%`)


# table 

tt1 = adjusted_growth %>%
  select(Country.Name2,argnp=est)

tt1 = prop_iota %>% 
  select(Country.Name2,propiota=est2) %>%
  right_join(tt1)

# iota
tt1 = summary(S_M6,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9,
         iota=paste0(sprintf('%.1f',`50%`*100),"% (",sprintf('%.1f',`2.5%`*100),"; ",sprintf('%.1f',`97.5%`*100),")")) %>%
  left_join(country_key) %>%
  select(Country.Name2,iota) %>%
  right_join(tt1)
  

# pred
tt1 = summary(S_M6,pars="F_output")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=rep(1:9,each=D_S_M6$L),
         year=rep(1999+(1:D_S_M6$L),9)) %>%
  left_join(country_key) %>%
  filter(year==2018) %>%
  mutate(CountryID=1:9,
         pred2018=paste0(sprintf('%.1f',`50%`*100),"% (",sprintf('%.1f',`2.5%`*100),"; ",sprintf('%.1f',`97.5%`*100),")")) %>%
  select(Country.Name2,pred2018) %>%
  right_join(tt1)

# tau
tt1 = summary(S_M6,pars="tau")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9,
         tau=paste0(signif(`50%`,2)," (",signif(`2.5%`,2),"; ",signif(`97.5%`,2),")")) %>%
  left_join(country_key) %>%
  select(Country.Name2,tau)%>%
  right_join(tt1)

# nu
tt1 = summary(S_M6,pars="nu")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9,
         nu=paste0(1999+signif(`50%`,1)," (",1999+signif(`2.5%`,1),"; ",1999+signif(`97.5%`,1),")")) %>%
  left_join(country_key) %>%
  select(Country.Name2,nu) %>%
  right_join(tt1)

# overall resistance
yy = rstan::extract(S_M6,pars="y")[[1]]
yy_2 = yy[,,19,2]
yy_3 = yy[,,19,3]
yy_4 = yy[,,19,4]
yy_5 = yy[,,19,5]
yy_6 = yy[,,19,6]

yy_res_overall = (yy_4 + yy_5 + yy_6)/(yy_2 + yy_3+yy_4+yy_5+yy_6)
tt1 = t(apply(yy_res_overall,2,qsum2)) %>%
  tbl_df() %>%
  mutate(CountryID=1:9,
         overallres=paste0(sprintf('%.1f',`50%`*100),"% (",sprintf('%.1f',`2.5%`*100),"; ",sprintf('%.1f',`97.5%`*100),")")) %>%
  left_join(country_key) %>%
  select(Country.Name2,overallres) %>%
  right_join(tt1)

select(tt1,Country=Country.Name2,shift=nu,intensity=tau,argnp,Prev2000=iota,Prev2018=pred2018,propiota,overallres2018=overallres) %>%
  arrange(Country)

select(tt1,Country=Country.Name2,shift=nu,intensity=tau,Prev2018=pred2018,Prev2000=iota,propiota,argnp) %>%
  arrange(Country) %>%
  xtable() %>%
  print(include.rownames=FALSE)

select(tt1,Country=Country.Name2,shift=nu,intensity=tau,Prev2018=pred2018,Prev2000=iota,propiota,argnp) %>%
  arrange(Country) %>%
  write.table(., file = paste0(fig_path,"table1.txt"), sep = ",", quote = FALSE, row.names = F)

# suppression

## Figure res among treated ~ PDR
yy = rstan::extract(S_M6,pars="y")[[1]]
yy_2 = yy[,,1:19,2]
yy_3 = yy[,,1:19,3]
yy_4 = yy[,,1:19,4]
yy_5 = yy[,,1:19,5]
yy_6 = yy[,,1:19,6]

yy_prop_supp = (yy_3+yy_6) / (yy_3+yy_5+yy_6)
yy_prop_supp_b = NULL
for(i in 1:19) {
  tmp = yy_prop_supp[,,i] %>%
    as.data.frame() %>%
    tbl_df() %>%
    gather("CountryID","supp",1:9) %>%
    mutate(CountryID=as.numeric(gsub("V","",CountryID)),
           year=i+1999) %>%
    group_by(CountryID) %>%
    mutate(it=row_number())
  yy_prop_supp_b = bind_rows(yy_prop_supp_b,tmp)
  cat(i)
}

yy_prop_supp_b %>%
  group_by(CountryID,year) %>%
  summarise(m=median(supp),l=quantile(supp,0.025),h=quantile(supp,0.975)) %>%
  mutate(qsum=paste0(round(m,3)," (",round(l,3),"; ",round(h,3),")")) %>%
  left_join(country_key) %>%
  filter(year==2018) %>%
  arrange(m)


## Figures ----------------------------------------------------------------------


## Figure indicators
plot_indicators(samples=S_M6,data.list=D_S_M6,key=country_key)
ggsave(file=paste0(fig_path,"fig2b.pdf"),width=8,height=5)


## Figure PDR
plot_pdr(samples=S_M6,data.list=D_S_M6,key=country_key,ylim = .3)
ggsave(file="figures/post_pdr_M6.pdf",width=5.5,height=4)


plot_pdr3(samples=S_M6,data.list=D_S_M6,key=country_key,ylim = .3,outliers=c(13,50))
ggsave(file="figures/fig4.pdf",width=5.5,height=4)

## Figure ART rollout
g3A = summary(S_M6,pars="tau_t")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  mutate(CountryID=rep(1:9,each=D_S_M6$L),
         year=rep(1999+(1:D_S_M6$L),9)) %>%
  left_join(country_key) %>%
  ggplot(.) +
  geom_line(aes(x=year,y=`50%`,colour=Country.Name2),size=.8) +
  labs(x="Year",y="Rate of ART initiation (per year)",colour=NULL) +
  # bwa 
  annotate(geom="segment",x=2002,y=.05,xend=2003.5,yend=.03) +
  annotate(geom="label",x=2002,y=.05,label="BWA",size=2.5) +
  # moz
  annotate(geom="segment",x=2012,y=.04,xend=2009.6,yend=.045) +
  annotate(geom="label",x=2012,y=.04,label="MOZ",size=2.5) +
  # nam
  annotate(geom="segment",x=2006.5,y=.16,xend=2008,yend=.14) +
  annotate(geom="label",x=2006.5,y=.16,label="NAM",size=2.5) +
  # lso
  annotate(geom="segment",x=2017,y=.07,xend=2016.5,yend=.09) +
  annotate(geom="label",x=2017,y=.07,label="LSO",size=2.5) +
  # zwe
  annotate(geom="segment",x=2010.5,y=.22,xend=2012.6,yend=.21) +
  annotate(geom="label",x=2010.5,y=.22,label="ZWE",size=2.5)  +
  # swz
  annotate(geom="segment",x=2012.5,y=.25,xend=2015.1,yend=.223) +
  annotate(geom="label",x=2012.5,y=.25,label="ESW",size=2.5)  +
  # mwi
  annotate(geom="segment",x=2017,y=.20,xend=2016.5,yend=.22) +
  annotate(geom="label",x=2017,y=.20,label="MWI",size=2.5)  +
  # zmb
  annotate(geom="segment",x=2017,y=.15,xend=2016.5,yend=.17) +
  annotate(geom="label",x=2017,y=.15,label="ZMB",size=2.5)  +
  # zaf
  annotate(geom="segment",x=2013.9,y=.11,xend=2012.2,yend=.122) +
  annotate(geom="label",x=2013.9,y=.11,label="RSA",size=2.5)  +
  theme(legend.position="none")
g3A

## Figure res among treated ~ PDR
yy = rstan::extract(S_M6,pars="y")[[1]]
yy_2 = yy[,,1:19,2]
yy_3 = yy[,,1:19,3]
yy_4 = yy[,,1:19,4]
yy_5 = yy[,,1:19,5]
yy_6 = yy[,,1:19,6]

yy_pdr=rstan::extract(S_M6,pars="F_output")[[1]]


yy_prop_res = (yy_4+yy_5+yy_6) / (yy_2+yy_3+yy_4+yy_5+yy_6)
yy_prop_res_b = NULL
for(i in 1:19) {
  tmp = yy_prop_res[,,i] %>%
    as.data.frame() %>%
    tbl_df() %>%
    gather("CountryID","res_among_res",1:9) %>%
    mutate(CountryID=as.numeric(gsub("V","",CountryID)),
           year=i+1999) %>%
    group_by(CountryID) %>%
    mutate(it=row_number())
  tmp2 = yy_pdr[,,i] %>%
    as.data.frame() %>%
    tbl_df() %>%
    gather("CountryID","pdr",1:9) %>%
    mutate(CountryID=as.numeric(gsub("V","",CountryID)),
           year=i+1999) %>%
    group_by(CountryID) %>%
    mutate(it=row_number())
  yy_prop_res_b = bind_rows(yy_prop_res_b,left_join(tmp,tmp2))
  cat(i)
}

yy_prop_res_c = yy_prop_res_b %>%
  group_by(CountryID,year) %>%
  summarise(res=median(res_among_res),pdr=median(pdr)) %>%
  left_join(country_key) %>%
  mutate(country2=paste0(Country.Name2," (",Country.Abb,")"))

yy_label = filter(yy_prop_res_c,CountryID==6) %>%
  mutate(x=pdr/10+.02,y=res/10+.25) 

g3B = ggplot() +
  geom_abline(slope=1,linetype=2,color="grey") +
  geom_point(data=yy_prop_res_c,aes(x=pdr,y=res,colour=country2)) +
  geom_line(data=yy_prop_res_c,aes(x=pdr,y=res,colour=country2)) +
  annotate("text",x=0.05,y=0.307,label="2018",size=2.5) +
  annotate("text",x=0.018,y=0.24,label="2000",size=2.5) +
  annotate("segment",x=0.021,y=0.28,xend=0.028,yend=0.295,arrow=arrow(length=unit(.1,"cm"))) +
  geom_point(data=yy_label,aes(x=x,y=y),size=.5) +
  geom_line(data=yy_label,aes(x=x,y=y)) +
  labs(x="NNRTI PDR",y="Overall NNRTI resistance",colour=NULL) +
  scale_x_continuous(labels=scales::percent,limits=c(0,.4))+
  scale_y_continuous(labels=scales::percent,limits=c(0,.4)) +
  theme(legend.position=c(.8,.3),
        legend.margin=margin(1,3,3,3,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(10,"pt"))
g3B
# plot_grid(g3A,g3B,g3C,g3D,ncol=2,labels=c("A","B","C","D"))

## Figure iomegas
kk = arrange(country_key,Country.Name2)
g3C = summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  ggplot(.) +
  geom_linerange(aes(x=Country.Name2,ymin=`25%`,ymax=`75%`,colour=Country.Name2),size=1.1) +
  geom_linerange(aes(x=Country.Name2,ymin=`2.5%`,ymax=`97.5%`,colour=Country.Name2),size=.5) +
  geom_point(aes(x=Country.Name2,y=`50%`,colour=Country.Name2),shape="-",size=7) +
  scale_x_discrete(labels=kk$Country.Abb) +
  labs(x="Country",y="TFNR") +
  coord_cartesian(ylim=c(0,.55)) +
  theme(legend.position="none")
g3C

# tt_omega_lims = group_by(tt,CountryID) %>%
#   mutate(low=quantile(omega,0.025),
#          high=quantile(omega,0.975)) %>%
#   filter(omega<high & omega>low)
# ggplot(tt_omega_lims) +
#   geom_violin(aes(x=Country.Name2,y=omega,colour=Country.Name2,fill=Country.Name2),weight=2) +
#   coord_cartesian(ylim=c(0,.55)) +
#   scale_x_discrete(labels=kk$Country.Abb) 
#   



# Map

oms = adjusted_growth %>%
  rename(region=Country.Name)

africa = map_data("world") %>% 
  tbl_df() %>%
  filter(region %in% c("Botswana","Zambia","Lesotho","Mozambique","Malawi","Zimbabwe","Namibia","Swaziland","South Africa",
                       "Angola","Democratic Republic of the Congo","Tanzania","Republic of Congo","Comores")) %>%
  left_join(oms)
meds = group_by(africa,Country.Abb) %>%
  summarise(xmed=median(long),ymed=median(lat)) %>%
  mutate(x2=xmed,y2=ymed) %>%
  filter(!is.na(Country.Abb))
meds[meds$Country.Abb=="LSO",2:3] = list(35,-32)
meds[meds$Country.Abb=="MOZ",2:5] = list(40,-21,36,-18)
meds[meds$Country.Abb=="MWI",2:5] = list(36,-7,34,-11)
meds[meds$Country.Abb=="ZAF",2:5] = list(24,-29,24,-29)
meds[meds$Country.Abb=="SWZ",2:5] = list(37,-27,31.5,-26.6)
meds[meds$Country.Abb=="ZMB",2:5] = list(26,-15,26,-15)
meds[meds$Country.Abb=="NAM",2:5] = list(17.1,-20.0,17.1,-20 )
meds[meds$Country.Abb=="BWA",2:5] = list(24.4 ,-22.3,24.4, -22.3)
meds = left_join(meds,country_key)

africa2 = filter(africa,Country.Abb=="LSO")

g3C = ggplot() +
  # theme(panel.background = element_rect(fill="lightblue")) +
  geom_polygon(data=africa,aes(x = long, y = lat, group = group, fill=`50%`), colour = "black",alpha=.7) +
  geom_polygon(data=africa2,aes(x = long, y = lat, group = group), fill="white", colour = "black",alpha=.7) +
  geom_polygon(data=africa2,aes(x = long, y = lat, group = group, fill=`50%`), colour = "black",alpha=.7) +
  geom_segment(data=meds,aes(x=xmed,y=ymed,xend=x2,yend=y2)) +
  geom_label(data=meds,aes(x=xmed,y=ymed,label=Country.Abb2),size=2.5) +
  coord_map(xlim=c(0,43),ylim=c(-36,-4)) +
  scale_fill_gradient(low="green",high="red",trans="log",na.value="white",breaks=c(0.001,0.01,0.1,1,10)) +
  labs(x="Longitude",y="Latitude",fill="Fragility index") +
  theme(legend.position=c(.15,.25),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(12,"pt"))
g3C



## Figure pred20108 ~ omega + cum

xx = group_by(tt,Country.Abb) %>%
  summarise(x=median(cum_art/cum_inc),y=median(argnt),col=median(pred2018)) %>%
  mutate(xlab=x+.04,
         ylab=y+.1)
xx[xx$Country.Abb=="MOZ",c("xlab","ylab")] = list(.28,.4)
xx[xx$Country.Abb=="LSO",c("xlab","ylab")] = list(.33,.8)
xx[xx$Country.Abb=="ZWE",c("xlab","ylab")] = list(.37,.6)
xx[xx$Country.Abb=="ZMB",c("xlab","ylab")] = list(.42,.5)
xx[xx$Country.Abb=="MWI",c("xlab","ylab")] = list(.5,.7)
xx[xx$Country.Abb=="BWA",c("xlab","ylab")] = list(.57,0)
xx[xx$Country.Abb=="NAM",c("xlab","ylab")] = list(.59,.5)
xx = left_join(xx,country_key)

g3D = ggplot() +
  geom_point(data=tt,aes(x=cum_art/cum_inc,y=argnt,colour=pred2018),alpha=.07)  + 
  geom_point(data=xx,aes(x=x,y=y),shape=1) +
  scale_color_gradient(low="green",high="red",trans="sqrt",labels=scales::percent) +
  # scale_y_log10() +
  scale_x_continuous(labels=scales::percent) +
  coord_cartesian(ylim=c(0,5)) +
  labs(x="Cumulative number of ART initiators",y="Fragility index",colour="NNRTI PDR\nin 2018") +
  geom_segment(data=xx,aes(x=x,y=y,xend=xlab,yend=ylab)) +
  geom_label(data=xx,aes(x=xlab,y=ylab,label=Country.Abb2),size=2.5) +
  theme(legend.position=c(.2,.72),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(12,"pt"))

g3D

plot_grid(g3A,g3C,g3D,ncol=1,labels=c("A","B","C","D"),rel_heights = c(1,1,1))
ggsave(file=paste0(fig_path,"fig3.pdf"),width=4,height=9.75)




tt$cum_art2 = tt$cum_art/tt$cum_inc
getmode = function(x) {
  dens = density(x)
  n = which(dens$y==max(dens$y))
  return(list(x=x[n],y=dens$y[n]))
}
plotpar = c("beta","tau","xi","nu","kappa","delta","omega","cum_art2","iota","pred2018")
plotparnames = c("Transmission rate","Max. treatment rate","Slope of ART R-O","Shift of ART R-O","Switching rate","Mortality rate",
                 "TFNR","Cum. ART initiators","Res. in 1999","Res. in 2018")
gglist = list()
for(i in 1:length(plotpar)) {
  for(j in 1:length(plotpar)) {
    if(i!=j) {
      xx = group_by(tt,Country.Abb) %>%
        summarise(x=median(get(plotpar[i])),y=median(get(plotpar[j]))) %>%
        mutate(CountryID2=c(1,3:7,2,8,9)) 
      gglist[[i+((j-1)*10)]] = ggplot() +
        geom_point(data=tt,aes_string(x=plotpar[i],y=plotpar[j]),color="gray50",alpha=.05)  + 
        labs(x=plotparnames[i],y=plotparnames[j]) 
    } else {
      xx = group_by(tt,Country.Abb) %>%
        summarise(x=getmode(get(plotpar[i]))[[1]],y=getmode(get(plotpar[i]))[[2]]) %>%
        mutate(CountryID2=c(1,3:7,2,8,9)) 
      yshift = max(xx$y)/15
      gglist[[i+((j-1)*10)]] = ggplot() +
        geom_density(data=tt,aes_string(x=plotpar[i],fill="Country.Name2"),alpha=.5)  + 
        scale_fill_discrete(guide=FALSE) +
        labs(x=plotparnames[i],y="Density",colour="NNRTI PDR\n(2018)")
    }
  }
}

# plot_grid(plotlist=gglist,ncol=9)
ggsave(plot=plot_grid(plotlist=gglist,ncol=10),file="figures/plot_cor_pars.pdf",width=15,height=15)


xx = group_by(tt,Country.Abb) %>%
  summarise(x=median(cum_inc),y=median(kappa),col=median(pred2018))
ggplot() +
  geom_point(data=tt,aes(x=cum_inc,y=kappa,colour=pred2018),alpha=.05)  + 
  geom_point(data=xx,aes(x=x,y=y),shape=1) +
  scale_color_gradient(low="green",high="red",trans="sqrt",labels=scales::percent) +
  labs(x="Cumulative number of ART initiators",y="Beta",colour="NNRTI PDR\n(2018)") +
  theme(legend.position=c(.88,.72),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(12,"pt"))

## Correlation


omegas = rstan::extract(S_M6,"omega")[[1]]
kappas = rstan::extract(S_M6,"kappa")[[1]]
argnp = omegas/kappas
country_char

cc = NULL
for(i in 1:4000) {
  dd = data.frame(y=argnp[i,],
                  CountryID=1:9,stringsAsFactors = FALSE) %>%
    tbl_df() %>%
    left_join(country_key) %>%
    left_join(country_char,by = "Country.Abb") %>%
    select(y,Indicator.Name,value) %>%
    spread(Indicator.Name,value)
  cc = rbind(cc,cor(dd,method="spearman")[1,-1])
  print(i)
}
t(apply(cc,2,qsum)) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(fin=paste0(round(`50%`,2)," (",round(`2.5%`,2),"; ",round(`97.5%`,2),")")) %>%
  arrange(`50%`) %>%
  select(rowname,fin) %>% 
  xtable() %>%
  print(include.rownames=FALSE)
  
t(apply(cc,2,qsum)) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(fin=paste0(round(`50%`,2)," (",round(`2.5%`,2),"; ",round(`97.5%`,2),")")) %>%
  arrange(`50%`) %>%
  select(rowname,fin)  %>%
  write.table(., file = paste0(fig_path,"table2.txt"), sep = ",", quote = FALSE, row.names = F)




  iotas = rstan::extract(S_M6,"iota")[[1]]
  country_char
  
  cc = NULL
  for(i in 1:4000) {
    dd = data.frame(y=iotas[i,],
                    CountryID=1:9,stringsAsFactors = FALSE) %>%
      tbl_df() %>%
      left_join(country_key) %>%
      left_join(country_char,by = "Country.Abb") %>%
      select(y,Indicator.Name,value) %>%
      spread(Indicator.Name,value)
    cc = rbind(cc,cor(dd,method="spearman")[1,-1])
    print(i)
  }
  t(apply(cc,2,qsum)) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(fin=paste0(round(`50%`,2)," (",round(`2.5%`,2),"; ",round(`97.5%`,2),")")) %>%
    arrange(`50%`) %>%
    select(rowname,fin) %>% 
    xtable() %>%
    print(include.rownames=FALSE)

summary(S_M6,pars="iota")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Abb,median_omega=`50%`,low_omega=`25%`,high_omega=`75%`) %>%
  right_join(country_char) %>%
  ggplot() +
  geom_pointrange(aes(x=norm_value,y=median_omega,ymax=high_omega,ymin=low_omega,colour=Country.Abb)) +
  facet_wrap(~Indicator.Name,scales="free")






omega = extract(S_M6,pars="omega")[[1]]
kappa = extract(S_M6,pars="kappa")[[1]]
mu = 1/(D_S_M6$life_expectancy-15)
beta = extract(S_M6,pars="beta")[[1]]
delta = extract(S_M6,pars="delta")[[1]]


Rn = omega
Rn2 = omega
for(i in 1:9) {
  Rn[,i] = omega[,i]/(omega[,i]+mu[i]) * (beta[,i]/(kappa[,i]+mu[i]+delta[,i]))
  Rn2[,i] = (1-exp(-omega[,i]*5)) * (beta[,i]/(kappa[,i]+mu[i]+delta[,i]))
}
t(apply(Rn2,2,qsum)) %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Name,everything( )) %>%
  arrange(mean) %>%
  ggplot() +
  geom_pointrange(aes(x=Country.Name2,y=`50%`,ymin=`2.5%`,ymax=`97.5%`)) 
