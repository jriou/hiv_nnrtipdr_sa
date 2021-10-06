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

l = load("models/M6_2019-11-18.Rdata")

check_hmc_diagnostics(S_M6)
ppM6 = c("beta","tau","nu","xi","eta","delta","kappa","mu_omega","sigma_omega","omega","mu_iota","sigma_iota","iota","sigma")
print(S_M6,pars=ppM6,digits_summary = 4)

print(S_M6,pars=c("mu_omega","sigma_omega","omega"),digits_summary = 4)
print(S_M6,pars=c("mu_iota","sigma_iota","iota"),digits_summary = 3)



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
  left_join(country_key)

# table 

tt1 = summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9,
         omega=paste0(signif(`50%`,2)," (",signif(`2.5%`,2),"; ",signif(`97.5%`,2),")"),
         prop_fail=paste0(sprintf('%.1f',(1-exp(-`50%`))*100),"% (",sprintf('%.1f',(1-exp(-`2.5%`))*100),"; ",sprintf('%.1f',(1-exp(-`97.5%`))*100),")")) %>%
  left_join(country_key) %>%
  select(Country.Name2,omega,prop_fail)


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

select(tt1,Country=Country.Name2,shift=nu,intensity=tau,TFNR=omega,Prev2000=iota,Prev2018=pred2018,overallres2018=overallres) %>%
  arrange(Country)

select(tt1,Country=Country.Name2,shift=nu,intensity=tau,TFNR=omega,Prev2000=iota,Prev2018=pred2018,overallres2018=overallres) %>%
  xtable() %>%
  print(include.rownames=FALSE)


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

data.list = D_S_M6
samples = S_M6
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

prop = list(A_prop = A_pred / D_pred *100,
            B_prop = B_pred / A_pred*100,
            C_prop = C_pred / A_pred*100,
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
  left_join(country_key) %>%
  bind_cols(obs=c(as.vector(t(data.list$A_est/data.list$D_est*100)),
                  as.vector(t(data.list$B_est/data.list$A_est*100)),
                  as.vector(t(data.list$C_est/data.list$A_est*100)),
                  as.vector(t(data.list$D_est/1e6)))) %>%
  filter(Country.Abb2=="RSA")

g_indic =
  ggplot(pred) +
  annotate("point",x=2000,y=0,colour=NA) +
  geom_point(aes(x=year,y=obs),colour="grey30",size=.4,shape=21,fill="white") +
  geom_ribbon(aes(x=year,ymin=`2.5%`,ymax=`97.5%`,fill=comp2),alpha=.2) +
  geom_line(aes(x=year,y=`50%`,colour=comp2),size=.6) +
  facet_wrap(~comp,
             nrow=1,
             scales = "free_y",
             labeller=labeller(comp=c("A_pred"="Prevalence (%)",
                                      "B_pred"="ART (%)",
                                      "C_pred"="Mortality (%)",
                                      "D_pred"="Adult pop. (million)"))) +
  scale_x_continuous(breaks=c(2000,2008,2016),labels=c("'00","'08","'16"),limits=c(1999,2018)) +
  scale_fill_brewer(palette="Spectral",guide=FALSE) +
  scale_colour_brewer(palette="Spectral",guide=FALSE) +
  labs(x="Year",y=NULL) 
g_indic

comp=c("F_output")
comp2 = c("NNRTI PDR (%)")
pred = summary(samples,pars="F_output")[[1]] %>%
  as.data.frame(.) %>%
  rownames_to_column() %>%
  tbl_df() %>%
  mutate(year=start+rep(1:L,G),
         CountryID=rep(1:G,each=L),
         comp=rep(comp,each=G*L),
         comp2=rep(comp2,each=G*L)) %>%
  left_join(country_key) %>%
  filter(Country.Abb2=="RSA")
cc = "purple"

obs = data.frame(k=data.list$F_k,
                 n=data.list$F_n,
                 CountryID=data.list$F_country,
                 year=data.list$F_t+start) %>%
  tbl_df() %>% 
  mutate(p=k/n,
         pinf=qbeta(0.025, k, n-k+1),
         psup=qbeta(0.975, k+1, n-k)) %>%
  left_join(country_key) %>%
  group_by(CountryID,year) %>%
  mutate(rank=row_number(),
         year2=year+(rank-1)*.2) %>%
  ungroup() %>%
  filter(Country.Abb2=="RSA")
g =
  ggplot() +
  geom_ribbon(data=pred,aes(x=year,ymin=`2.5%`*100,ymax=`97.5%`*100),fill=cc,alpha=.2) +
  geom_line(data=pred,aes(x=year,y=`50%`*100),alpha=.6,size=.4) +
  geom_pointrange(data=obs,aes(x=year2,y=p*100,ymin=pinf*100,ymax=psup*100),colour="grey30",size=.2,shape=21,fill="white",stroke=.7) +
  facet_wrap(~comp2) +
  scale_x_continuous(breaks=c(2002,2006,2010,2014,2018),labels=c("'02","'06","'10","'14","'18")) +
  labs(x="Year",y=NULL) 
g

plot_grid(g_indic,g,rel_widths=c(3,1),nrow=1)
ggsave(file="figures/poster_indicators_M6.pdf",width=7,height=2)


## Figure PDR
plot_pdr(samples=S_M6,data.list=D_S_M6,key=country_key,ylim = .3) +
  facet_wrap(~Country.Name2,ncol=5) 
ggsave(file="figures/poster_pdr_M6.pdf",width=7,height=3.2)


## Figure ART rollout
g3A = summary(S_M6,pars="tau_t")[[1]] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  mutate(CountryID=rep(1:9,each=D_S_M6$L),
         year=rep(1999+(1:D_S_M6$L),9)) %>%
  left_join(country_key) %>%
  ggplot(.) +
  geom_line(aes(x=year,y=`50%`,colour=Country.Name2),size=.8) +
  labs(x="Year",y="Rate of ART initiation (/y)",colour=NULL) +
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

# Map

oms = summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
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
meds[meds$Country.Abb=="LSO",2:3] = c(35,-32)
meds[meds$Country.Abb=="MOZ",2:5] = c(40,-21,36,-18)
meds[meds$Country.Abb=="MWI",2:5] = c(36,-7,34,-11)
meds[meds$Country.Abb=="ZAF",2:5] = c(24,-29,24,-29)
meds[meds$Country.Abb=="SWZ",2:5] = c(37,-27,31.5,-26.6)
meds[meds$Country.Abb=="ZMB",2:5] = c(26,-15,26,-15)
meds[meds$Country.Abb=="NAM",2:5] = c(17.1,-20.0,17.1,-20 )
meds[meds$Country.Abb=="BWA",2:5] = c(24.4 ,-22.3,24.4, -22.3)
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
  scale_fill_gradient(low="green",high="red",trans="log",na.value="white",breaks=c(0.001,0.01,0.1)) +
  labs(x="Longitude",y="Latitude",fill="TFNR") +
  theme(legend.position=c(.15,.28),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(12,"pt"))
g3C



## Figure pred20108 ~ omega + cum

xx = group_by(tt,Country.Abb) %>%
  summarise(x=median(cum_art/cum_inc),y=median(omega),col=median(pred2018)) %>%
  mutate(xlab=x+.04,
         ylab=y+.1)
xx[xx$Country.Abb=="MOZ",c("xlab","ylab")] = c(.28,.102)
xx[xx$Country.Abb=="LSO",c("xlab","ylab")] = c(.33,.12)
xx[xx$Country.Abb=="ZWE",c("xlab","ylab")] = c(.37,.11)
xx[xx$Country.Abb=="ZMB",c("xlab","ylab")] = c(.42,.085)
xx[xx$Country.Abb=="MWI",c("xlab","ylab")] = c(.46,.10)
xx[xx$Country.Abb=="BWA",c("xlab","ylab")] = c(.57,0)
xx[xx$Country.Abb=="NAM",c("xlab","ylab")] = c(.59,.05)
xx = left_join(xx,country_key)

g3D = ggplot() +
  geom_point(data=tt,aes(x=cum_art/cum_inc,y=omega,colour=pred2018),alpha=.07)  + 
  geom_point(data=xx,aes(x=x,y=y),shape=1) +
  scale_color_gradient(low="green",high="red",trans="sqrt",labels=scales::percent) +
  # scale_y_continuous(trans="sqrt") +
  scale_x_continuous(labels=scales::percent) +
  coord_cartesian(ylim=c(0,.55)) +
  labs(x="Cumulative number of ART initiators",y="Risk of NNRTI resistance (/y)",colour="NNRTI PDR\nin 2018") +
  geom_segment(data=xx,aes(x=x,y=y,xend=xlab,yend=ylab)) +
  geom_label(data=xx,aes(x=xlab,y=ylab,label=Country.Abb2),size=2.5) +
  theme(legend.position=c(.88,.8),
        legend.margin=margin(4,5,4,4,"pt"),
        legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
        legend.text=element_text(size=6.5),
        legend.key.height=unit(12,"pt"))
g3D

plot_grid(g3A,g3C,nrow=1,labels=c("A","B","C","D"),rel_heights = c(1,1,1))
ggsave(file="figures/poster_f3.pdf",width=7,height=2.8)

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
country_char

cc = NULL
for(i in 1:4000) {
  dd = data.frame(y=omegas[i,],
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


summary(S_M6,pars="omega")[[1]] %>%
  as.data.frame() %>% 
  tbl_df() %>%
  mutate(CountryID=1:9) %>%
  left_join(country_key) %>%
  select(Country.Abb,median_omega=`50%`,low_omega=`2.5%`,high_omega=`97.5%`) %>%
  right_join(country_char) %>%
  ggplot() +
  geom_pointrange(aes(x=norm_value,y=median_omega,ymax=high_omega,ymin=low_omega,colour=Country.Abb)) +
  facet_wrap(~Indicator.Name,scales="free")

