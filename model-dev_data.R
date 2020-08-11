
# Create merging table ---------------------------------------------------------------------------------------------------
country_key = data.frame(
  Country.Name = c("Botswana","Lesotho","Malawi","Mozambique","Namibia","South Africa","Swaziland","Zambia","Zimbabwe"),
  Country.Name2 = c("Botswana","Lesotho","Malawi","Mozambique","Namibia","South Africa","Eswatini","Zambia","Zimbabwe"),
  CountryID = 1:9,
  Country.Abb = c("BWA","LSO","MWI","MOZ","NAM","ZAF","SWZ","ZMB","ZWE"),
  Reg = c(1,1,1,1,1,1,1,1,1),
  Rollout = c(2002,2004,2003,2004,2004,2004,2003,2002,2004),
  life_expectancy = c(69,53,63,59,63,63,58,63,61),
  stringsAsFactors = FALSE
)


# PDR survey data from Gupta et al, 2017 -----------------------------------------------------------------------------------------------

gupta_raw = 
  read.csv("data/gupta_data.csv",sep=",",stringsAsFactors = FALSE,header=TRUE,dec=",") %>%
  as_tibble() %>%
  filter(Country != "") %>% 
  rename(Country.Name=Country) %>%
  filter(!is.na(n_nnrti),!is.na(n_geno))

# Filter and clean up Gupta data
gupta_data = gupta_raw %>%
  filter(Country.Name %in% country_key$Country.Name) %>%
  left_join(country_key) %>%
  arrange(CountryID) %>%
  mutate(pdr_nnrti=n_nnrti/n_geno,
         pdr_nnrti_low=qbeta(0.025, n_nnrti+.5, n_geno-n_nnrti+.5), # Bayesian credible interval using the Jeffreys prior using https://www.jstor.org/stable/2676784?seq=1#metadata_info_tab_contents
         pdr_nnrti_high=qbeta(0.975, n_nnrti+.5, n_geno-n_nnrti+.5),
         year=round(as.numeric(Mid_Point_Year))) %>%
  select(CountryID,Country.Name,Country.Name2,Country.Abb,Rollout,Study,dataset,year,n_geno,n_nnrti,pdr_nnrti,pdr_nnrti_low,pdr_nnrti_high) %>%
  group_by(CountryID,year) %>%
  mutate(yearshift=year+(row_number()-1)*.2) %>%
  left_join(country_key) 

# Updated data
updated_data = read.csv("litterature_update/extraction_julien.csv") %>%
  tbl_df() %>%
  filter(art_initiators_or_art_naive==1, adults==1, original_data==1) %>%
  left_join(country_key) %>%
  select(CountryID,Country.Name,Country.Name2,Country.Abb,Rollout,Study,dataset=source,year,n_geno,n_nnrti,life_expectancy) %>%
  mutate(pdr_nnrti=n_nnrti/n_geno,
         pdr_nnrti_low=qbeta(0.025, n_nnrti+.5, n_geno-n_nnrti+.5), # Bayesian credible interval using the Jeffreys prior using https://www.jstor.org/stable/2676784?seq=1#metadata_info_tab_contents
         pdr_nnrti_high=qbeta(0.975, n_nnrti+.5, n_geno-n_nnrti+.5),
         Reg=1) %>%
  bind_rows(gupta_data) %>%
  arrange(CountryID,year) %>%
  group_by(CountryID,year) %>%
  mutate(yearshift=year+(row_number()-1)*.2) %>%
  ungroup()

ggplot(updated_data) +
  geom_pointrange(aes(x=yearshift,y=pdr_nnrti,ymax=pdr_nnrti_high,ymin=pdr_nnrti_low)) +
  facet_wrap(~Country.Name2)


# Population data from World Bank -------------------------------------------------------------------------

wb_pop = read.csv("data/API_SP.POP.1564.TO_DS2_en_csv_v2_64923.csv",skip=4) %>%
  tbl_df() %>%
  left_join(country_key,by = c("Country.Name"="Country.Name2"))  %>%
  filter(!is.na(CountryID)) %>%
  gather("year","adult_pop",5:63) %>%
  mutate(year=as.numeric(gsub("X","",year))) %>%
  filter(year>=1999) %>%
  select(CountryID,year,adult_pop)

# ggplot(wb_pop) +
#   geom_point(aes(x=year,y=adult_pop)) +
#   facet_wrap(~Country.Name2)

# HIV prevalence from UNAIDS -----------------------------------------------------------------------------

unaids_prevalence = read.csv("data/People living with HIV_Number of people living with HIV_Population_ Adults (15+).csv",
         stringsAsFactors = FALSE,strip.white=TRUE,na.strings=c("...")) %>%
  tbl_df() %>%
  rename(Country.Name2=Country) %>%
  left_join(country_key)  %>%
  filter(!is.na(CountryID)) %>%
  select(-ends_with("upper"),-ends_with("lower")) %>%
  gather("year","adult_prevalence",2:30) %>%
  mutate(year=as.numeric(gsub("X","",year))) %>%
  filter(year>=1999) %>%
  arrange(CountryID,year)

# ggplot(unaids_prevalence) +
#   geom_point(aes(x=year,y=adult_prevalence)) +
#   facet_wrap(~Country.Name2)

# HIV treatment coverage from UNAIDS ---------------------------------------------------------------------------

unaids_art_coverage_2010_2018 = read.csv("data/Treatment cascade_People living with HIV receiving ART (#)_Population_ Adults (15+).csv",
         stringsAsFactors = FALSE,strip.white=TRUE,na.strings=c("...")) %>%
  tbl_df() %>%
  rename(Country.Name2=Country) %>%
  left_join(country_key)  %>%
  filter(!is.na(CountryID)) %>%
  gather("year","adult_art_coverage",2:10) %>%
  mutate(year=as.numeric(gsub("X","",year)))  %>%
  arrange(CountryID,year)

ggplot() +
  geom_point(data=unaids_art_coverage_2010_2018,aes(x=year,y=adult_art_coverage)) +
  facet_wrap(~Country.Name2)

## Complete with HNP: health and nutrition variables (World bank)

hnp = read.csv("data/country_covariables/HNP_StatsData.csv") %>%
  as.tibble() %>%
  select(1:4,35:62)
names(hnp)[5:32] = 1990:2017
hnp = gather(hnp,"year","value",5:32) %>%
  mutate(year=as.numeric(year)) %>%
  select(Country.Abb=Country.Code,year,Indicator.Code,Indicator.Name,value) %>%
  left_join(country_key) %>%
  filter(!is.na(CountryID),year<=2016) %>%
  arrange(CountryID,year) 
hnp_cov = filter(hnp,Indicator.Code %in% c("SH.HIV.ARTC.ZS")) %>%
  select(Country.Name,CountryID,year,adult_proportion_art2=value) %>%
  arrange(CountryID,year) %>%
  left_join(unaids_prevalence) %>%
  mutate(adult_art_coverage=adult_prevalence*adult_proportion_art2/100) %>%
  filter(year>=1999)

# ggplot() +
#   geom_point(data=unaids_art_coverage_2010_2018,aes(x=year,y=adult_art_coverage)) +
#   geom_point(data=hnp_cov,aes(x=year,y=adult_art_coverage),colour="red") +
#   facet_wrap(~Country.Name2)

# Combine
unaids_art_coverage = expand.grid(CountryID=country_key$CountryID,year=1999:2018) %>%
    tbl_df() %>%
    left_join(select(unaids_art_coverage_2010_2018,CountryID,year,adult_art_coverage)) %>%
    left_join(select(hnp_cov,CountryID,year,adult_art_coverage2=adult_art_coverage)) %>%
    mutate(adult_art_coverage=ifelse(is.na(adult_art_coverage),adult_art_coverage2,adult_art_coverage),
           adult_art_coverage=ifelse(is.na(adult_art_coverage),0,adult_art_coverage)) %>%
    select(CountryID,year,adult_art_coverage) %>%
    left_join(country_key) %>%
    arrange(CountryID,year) 
  
# ggplot() +
#   geom_point(data=unaids_art_coverage_2010_2018,aes(x=year,y=adult_art_coverage)) +
#   geom_point(data=hnp_cov,aes(x=year,y=adult_art_coverage),colour="red") +
#   geom_line(data=unaids_art_coverage,aes(x=year,y=adult_art_coverage)) +
#   facet_wrap(~Country.Name2)


# AIDS-related mortality from UNAIDS ---------------------------------------------------------------------------

unaids_mortality = read.csv("data/AIDS-related deaths_Number of AIDS-related deaths_Population_ Adults (15+).csv",
                             stringsAsFactors = FALSE,strip.white=TRUE,na.strings=c("...")) %>%
  tbl_df() %>%
  rename(Country.Name2=Country) %>%
  left_join(country_key)  %>%
  filter(!is.na(CountryID)) %>%
  select(-ends_with("upper"),-ends_with("lower")) %>%
  gather("year","adult_mortality",2:30) %>%
  mutate(year=as.numeric(gsub("X","",year))) %>%
  filter(year>=1999) %>%
  arrange(CountryID,year)

# ggplot(unaids_mortality) +
#   geom_point(aes(x=year,y=adult_mortality)) +
#   facet_wrap(~Country.Name2)

# Viral suppression data (from UNAIDS AIDSinfo spectrum) ---------------------------------------------------------------

unaids_cascade = read.csv("data/Treatment cascade_Progress towards 90-90-90 targets - All ages.csv",na.strings = c("..."," ...",""),
                              skip=1,stringsAsFactors = FALSE,strip.white=TRUE) %>%
  tbl_df() %>%
  mutate_all(function(x) gsub(">","",x)) %>%
  select(Country.Name2=Country,
         prop2015="Percent.of.people.on.ART.who.achieve.viral.suppression",
         prop2016="Percent.of.people.on.ART.who.achieve.viral.suppression.1",
         prop2017="Percent.of.people.on.ART.who.achieve.viral.suppression.2",
         prop2018="Percent.of.people.on.ART.who.achieve.viral.suppression.3") %>%
  left_join(country_key) %>%
  filter(!is.na(CountryID)) %>%
  gather("year","control_prop",2:5) %>%
  filter(!is.na(control_prop)) %>%
  mutate(year=as.numeric(gsub("prop","",year)),
         control_prop=as.numeric(control_prop)/100) 


# add suppression rate for mozambique from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5881256/#ihx021C6
unaids_cascade = bind_rows(unaids_cascade,data.frame(Country.Name="Mozambique",Country.Name2="Mozambique",CountryID=4,
                                                     Country.Abb="MOZ",Reg=1,Rollout=2004,year=2015,control_prop=.82))

ggplot(unaids_cascade) +
  geom_point(aes(x=year,y=control_prop,colour=Country.Abb)) +
  geom_line(aes(x=year,y=control_prop,colour=Country.Abb)) 

unaids_cascade_summary = arrange(unaids_cascade,CountryID,year) %>%
  group_by(CountryID) %>%
  summarise(control_prop=mean(control_prop))

# PMTCT data  (from UNAIDS) -----------------------------------------------------------------------------------------
unaids_pmtct = read.csv("data/Elimination of mother-to-child transmission_Pregnant women who received ARV for PMTCT.csv",
         na.strings = c("..."," ...",""),stringsAsFactors = FALSE,strip.white=TRUE) %>%
  tbl_df() %>%
  rename(Country.Name2=Country) %>%
  left_join(country_key) %>%
  filter(!is.na(CountryID)) %>%
  gather("year","pmtct",2:10) %>%
  mutate(year=as.numeric(gsub("X","",year))) %>%
  arrange(CountryID,year)

ggplot(unaids_pmtct) +
  geom_line(aes(x=year,y=pmtct,colour=Country.Abb)) +
  facet_wrap(~Country.Name2,scale="free")


# complete using PMTCT efficacy data  (from UNAIDS) -----------------------------------------------------------------------------------------

unaids_pmtct_efficacy = read_excel("data/HIV_estimates_from_1990-to-present.xlsx",
         skip=5,sheet = 1) %>%
  tbl_df() %>%
  rename(year=1,Country.Abb=2) %>%
  left_join(country_key) %>%
  filter(!is.na(CountryID)) %>%
  select(Country.Name,year,
         mothers_needing_PMTCT=`MothersNeedingAntiretroviralsForPreventingMother-to-ChildTransmission`,
         children_incidence=`Children(0-14)newlyinfectedwithHIV`) %>%
  mutate(mothers_needing_PMTCT=as.numeric(gsub("<| ","",mothers_needing_PMTCT)),
         children_incidence=as.numeric(gsub("<| ","",children_incidence)),
         ratio_infected=children_incidence/mothers_needing_PMTCT) %>%
  filter(year>=1999) %>%
  group_by(Country.Name) %>%
  mutate(scaled_ratio_infected=ratio_infected/max(ratio_infected),
         extrapolated_pmtct=mothers_needing_PMTCT*(1-scaled_ratio_infected))%>%
  left_join(country_key)

ggplot(unaids_pmtct_efficacy) +
  geom_line(data=unaids_pmtct_efficacy,aes(x=year,y=mothers_needing_PMTCT),colour="red") +
  geom_line(data=unaids_pmtct_efficacy,aes(x=year,y=children_incidence),colour="blue") +
  facet_wrap(~Country.Name,scale="free")

ggplot(unaids_pmtct_efficacy) +
  geom_line(data=unaids_pmtct_efficacy,aes(x=year,y=ratio_infected),colour="blue") +
  facet_wrap(~Country.Name,scale="free")

ggplot(unaids_pmtct_efficacy) +
  geom_line(data=unaids_pmtct_efficacy,aes(x=year,y=1-scaled_ratio_infected),colour="blue") +
  facet_wrap(~Country.Name,scale="free")

ggplot() +
  geom_line(data=unaids_pmtct_efficacy,aes(x=year,y=extrapolated_pmtct),colour="red") +
  geom_line(data=unaids_pmtct,aes(x=year,y=pmtct)) +
  facet_wrap(~Country.Name,scale="free")

# Combine
unaids_pmtct_extrapolated = expand.grid(CountryID=country_key$CountryID,year=1999:2018) %>%
  tbl_df() %>%
  left_join(select(unaids_pmtct,CountryID,year,pmtct)) %>%
  left_join(select(unaids_pmtct_efficacy,CountryID,year,extrapolated_pmtct)) %>% 
  mutate(pmtct=ifelse(is.na(pmtct),extrapolated_pmtct,pmtct)) %>%
  select(CountryID,year,pmtct) %>%
  left_join(country_key) %>%
  arrange(CountryID,year)



ggplot() +
  geom_point(data=unaids_pmtct_efficacy,aes(x=year,y=extrapolated_pmtct),colour="red") +
  geom_point(data=unaids_pmtct,aes(x=year,y=pmtct)) +
  geom_line(data=unaids_pmtct_extrapolated,aes(x=year,y=pmtct)) +
  facet_wrap(~Country.Name,scale="free_y")


# Manage country data -----------------------------------------------------------------------------------------------

## Merge everything
country_data = wb_pop %>%
  left_join(unaids_prevalence) %>%
  left_join(unaids_art_coverage) %>%
  left_join(unaids_mortality) %>%
  left_join(unaids_pmtct_extrapolated) %>%
  left_join(unaids_cascade)


# Extract country characteristics ----------- ------------------------------------------------------------------------------------

chosen_char = c("NY.GNP.PCAP.CD","SP.RUR.TOTL.ZS","SL.UEM.TOTL.ZS","SH.XPD.CHEX.PP.CD","SE.ADT.LITR.ZS","SH.MED.BEDS.ZS",
                "SH.XPD.OOPC.CH.ZS","SH.XPD.EHEX.CH.ZS")

tt_pmtct = filter(country_data,year %in% 2000:2018) %>%
  mutate(value=pmtct/adult_prevalence,
         Indicator.Code="PMTCT",Indicator.Name="PMTCT (% prevalence)") %>%
  select(Country.Abb,year,Indicator.Code,Indicator.Name,value) %>%
  left_join(country_key)

country_char = filter(hnp,Indicator.Code %in% chosen_char,year %in% 2000:2018)%>%
  bind_rows(tt_pmtct) %>%
  group_by(Country.Abb,Indicator.Code,Indicator.Name) %>%
  summarise(value=mean(value,na.rm=T)) %>%
  group_by(Indicator.Code) %>%
  mutate(mean=mean(value),
         sd=sd(value),
         norm_value=(value-mean)/sd)

ggplot(country_char) +
  geom_col(aes(x=Country.Abb,y=norm_value)) +
  facet_wrap(~ Indicator.Name,scales="free",ncol=3)

country_char_matrix = matrix(country_char$norm_value,ncol=length(chosen_char)+1,byrow = TRUE)
