##############################################################################
# Name of file: 02_Modelling TND
# Original author(s): Thiago Cerqueira
# Original date: 11 Nov 21
# Latest update author (if not using version control) - thiago.c.silva@fiocruz.br
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 4.1.1
# Description of content: Modelling
# Approximate run time: Unknown
##############################################################################
setwd("/dados/Analysis/Thiago/Scripts/CoronaVac-Waning/TND")
pacman::p_load(tidyverse,twang,tidylog,tidytable,tableone,biglm,mgcv,parallel,ipw,WeightIt,cobalt,gtsummary,rms, tictoc, MatchIt, survival)


vac_test_clean <- read_rds("vac_test_clean_0812.RDS")

# Datasets infection
vac_test <- db_tnd_analysis %>% mutate(idade_10=cut(idade, breaks=c(18,seq(30,80,by=10),Inf), right = FALSE,
                                                                                           labels=c(paste(c(18,seq(30,70, by=10)),seq(29,79, by=10),sep="-"),"80+"), ordered_result=T,
                                                                                           include.lowest=T),
                                                                      idade_cat=cut(idade, breaks=c(18,60,80,+Inf), right = FALSE,
                                                                                   labels=c(paste(c(18,60),c(59,79),sep="-"),"80+"), ordered_result=T,
                                                                                   include.lowest=T))


vac_test_rt_pcr <- vac_test %>% filter(tipo_teste=="RT-PCR")

vac_test_rt_pcr_severe <- vac_test%>% filter(tipo_teste=="RT-PCR") %>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()

# Define functions
models_fit <- function(dataset){
  if(min(dataset$idade)<56)
    bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                                    family = binomial, nthreads = 24, discrete = T,data = dataset)
  else
      bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
              family = binomial, nthreads = 24, discrete = T,data = dataset)
}
extract_ve <- function(model, name, booster) {
  round(exp(cbind("Odds Ratio" = coef(model), confint.default(model, level = 0.95))), digits = 3) %>%
    as.data.frame() %>%
    rownames_to_column(var = "term") %>%
    filter(str_detect(term, name)) %>%
    filter(case_when(str_detect(term,"_v3_")~str_detect(term,booster),
                     TRUE~TRUE)) %>% 
    mutate(term = str_remove(term, paste0(name,"_"))) %>%
    mutate(across(c(2:4), ~ (1 - .x) * 100, .names = "{.col}_ve")) %>%
    select(1, 5:7) %>%
    mutate(ve_ci = paste(sprintf("%.1f", `Odds Ratio_ve`), "% (", sprintf("%.1f",`97.5 %_ve`), "-", sprintf("%.1f",`2.5 %_ve`), ")")) %>% 
    mutate(term=factor(term,levels=c("v1_0:6","v1_7:13","v1_14+","v2_0:13",
                                     "v2_14:30","v2_31:60","v2_61:90","v2_91:120",
                                     "v2_121:150","v2_151:180","v2_181+","v3_0:6_BNT162b2",
                                     "v3_7:13_BNT162b2","v3_14:30_BNT162b2","v3_31+_BNT162b2"))) %>% arrange(term) 
}


vac_df_10_infection <- vac_test_rt_pcr %>% nest(-idade_10) %>%  mutate(fit = map(data,models_fit ))%>% mutate(ve=map(fit,extract_ve)) 
vac_df_cat_infection <- vac_test_rt_pcr %>% nest(-idade_cat) %>%  mutate(fit = map(data,models_fit ))%>% mutate(ve=map(fit,extract_ve)) 


my_tidy <- function(x, exponentiate =  TRUE, conf.level = 0.95, ...) {
  dplyr::bind_cols(
    broom::tidy(x, exponentiate = exponentiate, conf.int = FALSE),
    # calculate the confidence intervals, and save them in a tibble
    stats::confint.default(x, level = 0.95) %>%
      tibble::as_tibble() %>%
      rlang::set_names(c("conf.low", "conf.high"))  )
}

vac_df_10_severe <- vac_test_rt_pcr_severe %>% nest(-idade_10) %>%  mutate(fit = map(data,models_fit ))%>% mutate(ve=map(fit,extract_ve)) 
vac_df_cat_severe <- vac_test_rt_pcr_severe %>% nest(-idade_cat) %>%  mutate(fit = map(data,models_fit ))

vac_df_cat_severe%>% mutate(ve=map(fit,extract_ve,name="vs_type2CV",booster="BNT"))  %>%  unnest(ve) %>% select(-c(data,fit)) %>% View()


teste <- vac_df_cat_severe%>% mutate(events=map(data ,.%>% select(confirmado,vs_type2) %>% filter(str_detect(vs_type2,"CV_|uv"))  %>%  droplevels() %>%  tbl_summary(by=confirmado)%>%add_overall(last = T)))  %>%  select(events) 


vac_test_symptomatic <- vac_test_clean %>%filter(tipo_teste=="RT-PCR", assintomatico==0)

my_tidy <- function(x, exponentiate =  FALSE, conf.level = 0.95, ...) {
  dplyr::bind_cols(
    broom::tidy(x, exponentiate = exponentiate, conf.int = FALSE),
    # calculate the confidence intervals, and save them in a tibble
    stats::confint.default(x) %>%
      tibble::as_tibble() %>%
      rlang::set_names(c("conf.low", "conf.high"))  )
}





# Datasets severe
vac_test_severe <- vac_test_clean %>% filter(outcome_confirmado!="Outpatient_Positive", tipo_teste=="RT-PCR") %>% droplevels()
vac_test_severe_18 <- vac_test_severe %>% filter(idade<60)
vac_test_severe_60_80 <- vac_test_severe %>% filter(idade>59,idade<80)
vac_test_severe_80 <- vac_test_severe %>% filter(idade>79)


# Sensitivity 
vac_test_infection_sensi <- vac_test_clean 
vac_test_severe_sensi <- vac_test_clean %>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()



# Infection



map(list(vac_test_18_39,vac_test_40_59),models_up60)

m1_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                    family = binomial, data = vac_test_infection, nthreads = 24, discrete = T)
m18_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                     family = binomial, data = vac_test_infection_18, nthreads = 16, discrete = T)
m18_39_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                     family = binomial, data = vac_test_18_39, nthreads = 16, discrete = T)
m40_59_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                     family = binomial, data = vac_test_40_59, nthreads = 16, discrete = T)
m60_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                     family = binomial, data = vac_test_infection_6_8, nthreads = 16, discrete = T)
m80_infection <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                     family = binomial, data = vac_test_infection_80, nthreads = 16, discrete = T)

m1_symp <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                    family = binomial, data = vac_test_symptomatic, nthreads = 24, discrete = T)
m18_symp<- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                     family = binomial, data = vac_test_symp_18, nthreads = 16, discrete = T)
m60_symp <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                     family = binomial, data = vac_test_symp_6_8, nthreads = 16, discrete = T)
m80_symp <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                     family = binomial, data = vac_test_symp_80, nthreads = 16, discrete = T)







m1_severe <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                       family = binomial, data = vac_test_severe, nthreads = 16, discrete = T)
m1_severe_18 <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                    family = binomial, data = vac_test_severe_18, nthreads = 16, discrete = T)

m1_severe_6_8 <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                    family = binomial, data = vac_test_severe_60_80, nthreads = 16, discrete = T)

m1_severe_80 <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
                         family = binomial, data = vac_test_severe_80, nthreads = 16, discrete = T)


m1_infection_sensitivity <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                                family = binomial, data = vac_test_infection_sensi, nthreads = 16, discrete = T)
m1_severe_sensitivity <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                             family = binomial, data = vac_test_severe_sensi, nthreads = 16, discrete = T)

lapply(ls(pattern = "m"), function(x) {
  save(list=x, file = paste0(x, ".Rda"))
})

saveRDS(m1_severe_sensitivity,"Results/m1_seve_sensi.RDS")
## Tables
t1 <- gtsummary::tbl_regression(m1_severe, exponentiate=T,include=vs_type2, method=gam) %>% add_n(location="level") %>% add_nevent(location="level") %>% bold_labels() %>% modify_caption(caption = "Hospitalization or Death")
tidy_gam(m1_severe, exponentiate = T, conf.int = T) %>% filter(str_detect(term,"vs_type2CV")) %>%  View()
saveRDS(m1_severe_80,"Results/m1_severe_80.RDS")


teste <- round(exp(cbind("Odds Ratio" = coef(m1_symp), confint.default(m1_symp, level = 0.95))), digits = 3)%>%
  as.data.frame() %>% rownames_to_column(var="term") %>%
  filter(str_detect(term,"vs_type2"))  %>% 
  mutate(across(c(2:4), ~ (1 - .x)*100, .names = "{.col}_ve")) %>% select(1,5:7) 

extract_ve(m1_severe_80)


models_infection <- list(m1_infection,m1_infection_sensitivity,m1_symp, m18_infection,m60_infection,m80_infection)
names(models_infection) <- c("Infection","Sensitivity", "Symptomatic","18-59","60-79","80+")
map_dfr(models_infection,extract_ve,.id="Model") %>% 
  mutate(term=factor(term,levels=c("v1_0:6","v1_7:13","v1_14+","v2_0:13",
                                   "v2_14:30","v2_31:60","v2_61:90","v2_91:120",
                                   "v2_121:150","v2_151:180","v2_181+","v3_0:7_BNT162b2",
                                   "v3_8:13_BNT162b2","v3_14:30_BNT162b2","v3_31+_BNT162b2"))) %>% arrange(Model,term) %>% View()


models_severe <- list(m1_severe,m1_severe_sensitivity,m1_severe_18, m1_severe_6_8,m1_severe_80)
names(models_severe) <- c("Hosp_Death","Hosp_Death_Sensi", "18-59","60-79","80+")
map_dfr(models_severe,extract_ve,.id="Model") %>% 
  mutate(term=factor(term,levels=c("v1_0:6","v1_7:13","v1_14+","v2_0:13",
                                   "v2_14:30","v2_31:60","v2_61:90","v2_91:120",
                                   "v2_121:150","v2_151:180","v2_181+","v3_0:7_BNT162b2",
                                   "v3_8:13_BNT162b2","v3_14:30_BNT162b2","v3_31+_BNT162b2"))) %>% arrange(Model,term) %>% View()

### Sensitivity

dt_match <- vac_test_infection %>% mutate(id_rand=row_number())
dt_match_minimal <- dt_match %>% select(id_rand,confirmado,sexo,cod_ibge,idade,dt_coleta)
  
  
  
max_value <- 0
data_matched <- data.frame()
for (i in unique(dt_match_minimal$cod_ibge)){
  dta_m <- data.frame()
  dataset <- dt_match_minimal %>% filter(cod_ibge==i)
  try(m.out3 <- matchit(confirmado=="Confirmado" ~ idade + sexo  + dt_coleta, 
                    data=dataset,
                    caliper = c(dt_coleta = 3,idade=5),
                    std.caliper = c(FALSE,FALSE),
                    exact=c("sexo"),
                    distance = "glm",
                    replace=T,
                    ratio = 4))
  try(dta_m <- get_matches(m.out3))
  dta_m <- dta_m %>% mutate(subclass=as.numeric(subclass),
                            subclass=subclass+max_value)
  if (exists("data_matched")) data_matched <- bind_rows(data_matched,dta_m) else data_matched <- dta_m
  max_value <- max(as.numeric(data_matched$subclass))
}

data_matched2 <- data_matched %>% left_join(dt_match %>% select(-c(confirmado,sexo,cod_ibge,idade,dt_coleta)),by="id_rand")
saveRDS(data_matched2,"matched_.rds")


sensitivity_cases_h_d <- matched_  %>% 
  filter(outcome=="Hosp_Death" & confirmado=="Confirmado") %>% 
  dplyr::select(subclass)

data_matched_sensitivity <-
  matched_  %>% 
  mutate(h_d_pairs = case_when(
    subclass %in% sensitivity_cases_h_d$subclass ~ 1))

teste2 <- clogit(confirmado=="Confirmado"~ vs_type2+ prev_infected+n_comorb+ puerpera + gestante+strata(subclass), data=data_matched_sensitivity)

gtsummary::tbl_regression(teste,exponentiate=T)
### Descriptives

vac_test_infection %>%
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  mutate(
    age_group = case_when(
      idade < 60 ~ "18-59",
      idade < 80 ~ "60-79",
      idade > 79 ~ "80+"
    ),
    vs_type3 = case_when(str_detect(vs_type2, "CV_")~ as.character(vs_type2),
                         vs_type2=="uv"~"Unvaccinated",
                         TRUE~"Other Vaccines")
  ) %>% mutate(vs_type3=fct_relevel(factor(vs_type3),"Unvaccinated","CV_v1_0:6",
                                    "CV_v1_7:13","CV_v1_14+","CV_v2_0:13",
                                    "CV_v2_14:30","CV_v2_31:60","CV_v2_61:90",
                                    "CV_v2_91:120","CV_v2_121:150","CV_v2_151:180",
                                    "CV_v2_180+")) %>% 
  select(age_group, vs_type3, confirmado) %>%
  tbl_strata(
    strata = age_group,
    .tbl_fun =
      ~ .x %>%
        tbl_summary(by = confirmado) %>%
        add_n()
  )%>% bold_labels() 


vac_test_infection  %>%
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  mutate(
    age_group = case_when(
      idade < 60 ~ "18-59",
      idade < 80 ~ "60-79",
      idade > 79 ~ "80+"
    ),
    vs_type3 = case_when(
      str_detect(vs_type2, "CV_") ~ as.character(vs_type2),
      vs_type2 == "uv" ~ "Unvaccinated",
      TRUE ~ "Other Vaccines"
    ),
    vs_type3=case_when(
      vs_type3=="CV_v3_CV"~"Other Vaccines",
      vs_type3=="CV_v3_Ad26"~"Other Vaccines",
      vs_type3=="CV_v3_AZ"~"Other Vaccines",
      TRUE~vs_type3
    )
  ) %>%
  mutate(
    macro_region = str_extract(cod_ibge, "^.{1}"),
    macro_region = case_when(
      macro_region == "1" ~ "North",
      macro_region == "2" ~ "Northeast",
      macro_region == "3" ~ "Southeast",
      macro_region == "4" ~ "Central-west",
      macro_region == "5" ~ "South"
    )
  ) %>%
  mutate(vs_type3 = fct_relevel(
    factor(vs_type3), "Unvaccinated", "CV_v1_0:6",
    "CV_v1_7:13", "CV_v1_14+", "CV_v2_0:13",
    "CV_v2_14:30", "CV_v2_31:60", "CV_v2_61:90",
    "CV_v2_91:120", "CV_v2_121:150", "CV_v2_151:180",
    "CV_v2_180+"
  ),
  cv_time=case_when(
    str_detect(vs_type3,"CV_v2")~"Second dose-CV",
    str_detect(vs_type3,"CV_v1")~"First dose-CV",
    str_detect(vs_type2,"CV_v3_BNT162b2")~"Booster-CV BNT162"
  )) %>%
  select(idade, sexo,cv_time, raca,unique_id, age_group,tipo_teste, macro_region, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, n_comorb, prev_infected, hosp_event, death_event,outcome,tipo_teste, vs_type3, confirmado) %>%
  tbl_summary(by = confirmado,digits = list(all_categorical() ~ c(0, 1)),value = list(gestante:drc ~ "1", unique_id~"1"))%>% 
  bold_labels() %>% add_overall(last=T) %>% add_n()




vac_test_infection%>%  filter(str_detect(vs_type2,"CV_v3")) %>% droplevels() %>% 
  mutate(
    age_group = case_when(
      idade < 60 ~ "18-59",
      idade < 80 ~ "60-79",
      idade > 79 ~ "80+"
    ),
    vs_type3 = case_when(
      str_detect(vs_type2, "CV_") ~ as.character(vs_type2),
      vs_type2 == "uv" ~ "Unvaccinated",
      TRUE ~ "Other Vaccines"
    ),
    vs_type3=case_when(
      vs_type3=="CV_v3_CV"~"Other Vaccines",
      vs_type3=="CV_v3_Ad26"~"Other Vaccines",
      vs_type3=="CV_v3_AZ"~"Other Vaccines",
      TRUE~vs_type3
    )
  ) %>%
  mutate(
    macro_region = str_extract(cod_ibge, "^.{1}"),
    macro_region = case_when(
      macro_region == "1" ~ "North",
      macro_region == "2" ~ "Northeast",
      macro_region == "3" ~ "Southeast",
      macro_region == "4" ~ "Central-west",
      macro_region == "5" ~ "South"
    )
  )  %>% filter(vs_type2=="CV_v3_BNT162b2") %>% 
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  select(age_group,vs_type3, length,idade, sexo, raca,unique_id, age_group,tipo_teste, macro_region, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, n_comorb, prev_infected, hosp_event, death_event,outcome,tipo_teste, vs_type3, confirmado) %>%
      tbl_summary(by = vs_type3,digits = list(all_categorical() ~ c(0, 1)))%>% add_overall(last=T)%>% bold_labels() 




# Severe
vac_test_severe_sensitivity <- vac_test_model %>% mutate(idade=if_else(idade>99,100L,idade))%>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()
vac_test_severe <- vac_test_model %>% filter(outcome_confirmado!="Outpatient_Positive", tipo_teste=="RT-PCR") %>% droplevels()



#Only CV-BNT162
m1_v3 <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+ prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
             family = binomial, data = vac_test_infection, nthreads = 16)
m1_v3_sensi <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+ prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
             family = binomial, data = vac_test_model, nthreads = 16)
saveRDS(m1,"m1_infection.RDS")


#CV-Bnt and CV-CV
m1_v3_sensi_two <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+ prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
             family = binomial, data = vac_test_infection_sensi, nthreads = 24)


m1 <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+ prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
          family = binomial, data = vac_test_infection, nthreads = 16)



m1_sensi <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+ prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
          family = binomial, data = vac_test_sensitivity, nthreads = 10)
saveRDS(m1_sensi,"m1_infection_sensi.RDS")



m2 <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
          family = binomial, data = vac_test_severe, nthreads = 10)
saveRDS(m2,"m2_severe.RDS")
rm(m2)
gc()
m2_sensi <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
          family = binomial, data = vac_test_severe_sensitivity, nthreads = 10)
saveRDS(m2_sensi,"m2_severe_sensi.RDS")
rm(m2_sensi)
gc()
m3 <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + puerpera + gestante + drc,
                family = binomial, data = vac_test_severe_18, nthreads = 10)
saveRDS(m3,"severe_18plus.RDS")
rm(m3)
gc()
m4 <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca+drc,
                family = binomial, data = vac_test_severe_60, nthreads = 10)
saveRDS(m4,"severe_60plus.RDS")
m5 <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + drc,
          family = binomial, data = vac_test_severe_60_80, nthreads = 10)
saveRDS(m5,"severe_60_79.RDS")
m6 <- bam(outcome_confirmado=="Hosp_Death_Positive" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+prev_infected+ sexo +uf + diabetes + obesidade + imunossupressao + dcardiaca + drc,
                family = binomial, data = vac_test_severe_80, nthreads = 10)
saveRDS(m6,"severe_80plus.RDS")





vac_test_model <- vac_test_model %>% mutate(vaccinated_unvacc=if_else(vs_type=="uv","uv","v"),
                                          cv_other=case_when(
                                            str_detect(vs_type,"CV")~"CoronaVac",
                                            vs_type!="uv"~ "other vaccines"
                                          ))

vac_test_infection <- vac_test_sensitivity %>% mutate(vaccinated_unvacc=if_else(vs_type=="uv","uv","v"),
                                            cv_other=case_when(
                                              str_detect(vs_type,"CV")~"CoronaVac",
                                              vs_type!="uv"~ "other vaccines",
                                              vs_type=="uv"~"Unvaccinated"
                                            ),
                                            cv_time=case_when(
                                              str_detect(vs_type,"CV")~as.character(vs_type),
                                              !str_detect(vs_type,"CV") & vs_type!="uv" ~ "other vaccines",
                                              vs_type=="uv"~"Unvaccinated"
                                            ))

vac_test_clean %>%
  select(idade, sexo, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, raca, n_comorb, prev_infected, tipo_teste, confirmado,hosp_event,death_event,outcome) %>%
  tbl_summary(by = confirmado) %>%
  bold_labels()


vac_test_infection %>%
  select(cv_time,confirmado) %>%
  tbl_summary(by = confirmado) %>%
  bold_labels() %>% add_overall(last=T)




