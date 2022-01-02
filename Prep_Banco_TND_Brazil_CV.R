library(tidyverse)
library(sparklyr)
spark_disconnect(sc)
unixtools::set.tempdir("/dados/tmp/") #Set temp dir (H2O)
conf <- sparklyr::spark_config()
# def temp dir_spark
spark_dir = "/dados/tmp/"
conf$`sparklyr.shell.driver-java-options` <-  paste0("-Djava.io.tmpdir=", spark_dir)
conf$`sparklyr.shell.driver-memory` <- "222G" # Memória
conf$spark.driver.maxResultSize <- "8g" #Memória - table
conf["sparklyr.connect.cores.local"] <- 32  #Número de cores
sc <- spark_connect(master = "local", version = "3.1.1",config = conf)


banco <- spark_read_parquet(sc,
                            "/dados/vac_srag_sg_11112021_v6.parquet/", 
                            memory = F)

# lista_ibp <- spark_read_csv(sc,"/dados/BDI_Municipalities-Level_Short.csv")
# lista_ibp <- lista_ibp %>% mutate(codmun6=as.character(codmun6))
# banco <- banco %>%  left_join(lista_ibp,by=c("d_codibge"="codmun6"))
############### TND
banco_prep <- banco %>%
  mutate( # criar uma variável com apenas astrazeneca em cada dose
    d1_nome_vacina = if_else(d1_nome_vacina %in% c("Covid-19-AstraZeneca", "Vacina Covid-19 - Covishield"), "Covid-19-AstraZeneca", d1_nome_vacina),
    d2_nome_vacina = if_else(d2_nome_vacina %in% c("Covid-19-AstraZeneca", "Vacina Covid-19 - Covishield"), "Covid-19-AstraZeneca", d2_nome_vacina))

tbl_vars(banco_prep)

bancosg <- banco_prep %>% select(id_vigvac,n1_sg_dt_datanotificacao:n5_sg_dt_dataencerramento) %>% filter(!is.na(n1_sg_dt_datanotificacao)|
                                                                                                            !is.na(n2_sg_dt_datanotificacao)|
                                                                                                            !is.na(n3_sg_dt_datanotificacao)|
                                                                                                            !is.na(n4_sg_dt_datanotificacao)|
                                                                                                            !is.na(n5_sg_dt_datanotificacao))
bancosrag <- banco_prep %>% select(id_vigvac,n1_srag_dt_notif:n2_srag_idade_srag) %>% filter(!is.na(n1_srag_dt_notif)|
                                                                                               !is.na(n2_srag_dt_notif))
# vacina
bancovac <- banco_prep %>%
  select(id_vigvac, d_sexo:d1_data_aplicacao, d1_vacina_categoria_nome, d2_nome_vacina:d2_data_aplicacao, d3_nome_vacina:d3_data_aplicacao) %>%
  filter(!is.na(d1_data_aplicacao))
spark_write_parquet(bancovac,"/dados/Analysis/Thiago/DB_processed/db_vacine111121/")


###
srag_long <- bancosrag %>%pivot_longer(cols = -id_vigvac, 
                                       names_to = c("level",".value"),
                                       values_to = "values",
                                       names_pattern = "(n[1|2])_(.*)")
srag_long <- srag_long %>% mutate(level=regexp_replace(level, "n", "")) %>% filter(!is.na(srag_dt_notif)) #remover letra n da variável level.


sg_long <- bancosg %>%pivot_longer(cols = -id_vigvac, 
                                   names_to = c("level",".value"),
                                   values_to = "values",
                                   names_pattern = "(n[1-5])_(.*)")

sg_long <- sg_long %>% mutate(level=regexp_replace(level, "n", "")) %>% filter(!is.na(sg_dt_datanotificacao))

sg_flowchart <- sg_long %>%
  select(id_vigvac, sg_dt_datanotificacao) %>%
  rename(dt_notificacao = sg_dt_datanotificacao) %>%
  mutate(sistema = "sg")
srag_flowchart <- srag_long %>%
  select(id_vigvac, srag_dt_notif) %>%
  rename(dt_notificacao = srag_dt_notif) %>%
  mutate(sistema = "srag")
srag_sg <- sdf_bind_rows(sg_flowchart,srag_flowchart)
spark_write_parquet(srag_sg,"/dados/Analysis/Thiago/DB_processed/flowchart_notifications/")

###Extract data about comorbidities and recode name to posterior bind with SRAG/ Select only necessary columns
sg_long <- sg_long %>% 
  mutate(sg_confirmado=case_when(sg_tipoteste_2 %in% c("Antígeno","RT-PCR") & 
                                   sg_resultadoteste %in% c("Positivo","Detectável","Reagente")~"Confirmado",
                                 sg_tipoteste_2 %in% c("Antígeno","RT-PCR") & 
                                   sg_resultadoteste %in% c("Negativo","Não detectável","Não Reagente") & 
                                   datediff(sg_dt_datateste,sg_dt_datainiciosintomas)<11~"Negativo")) %>% 
  ungroup() %>%
  filter(!is.na(sg_confirmado)) %>%
  mutate(
    sg_condicoes = if_else(is.na(sg_condicoes), "Vazio", sg_condicoes),
    assintomatico = if_else(instr(sg_sintomas, "Assintomático") > 0, 1L, 0L),
    diabetes = if_else(instr(sg_condicoes, "Diabetes") > 0, 1L, 0L),
    obesidade = if_else(instr(sg_condicoes, "Obesidade") > 0, 1L, 0L),
    imunossupressao = if_else(instr(sg_condicoes, "Imunossupressão") > 0, 1L, 0L),
    dcardiaca = if_else(instr(sg_condicoes, "cardíacas") > 0, 1L, 0L),
    puerpera = if_else(instr(sg_condicoes, "parto") > 0, 1L, 0L),
    drc = if_else(instr(sg_condicoes, "renais") > 0, 1L, 0L),
    gestante = if_else(instr(sg_condicoes, "Gestante") > 0, 1L, 0L),
    sg_racacor_2 = case_when(
      sg_racacor_2 == "Branca" ~ 1L,
      sg_racacor_2 == "Preta" ~ 2L,
      sg_racacor_2 == "Amarela" ~ 3L,
      sg_racacor_2 == "Parda" ~ 4L,
      sg_racacor_2 == "Indigena" ~ 5L,
      is.na(sg_racacor_2) ~ NA_integer_
    ),
    sg_sexo_2 = case_when(
      sg_sexo_2 == "Indefinido" ~ "I",
      sg_sexo_2 == "Feminino" ~ "F",
      sg_sexo_2 == "Masculino" ~ "M",
      is.na(sg_sexo_2) ~ NA_character_
    )
  ) %>%
  rename(
    dt_notificacao = sg_dt_datanotificacao,
    dt_inicio_sintomas = sg_dt_datainiciosintomas,
    dt_coleta = sg_dt_datateste,
    sexo = sg_sexo_2,
    raca = sg_racacor_2,
    idade = sg_idade,
    confirmado = sg_confirmado,
    cod_ibge = sg_codigo_ibge_resid,
    tipo_teste = sg_tipoteste_2
  ) %>%
  mutate(sistema = "sg") %>%
  # adicionar variável indicando o tipo de notificação
  select(
    id_vigvac,
    dt_notificacao,
    dt_inicio_sintomas,
    dt_coleta,
    idade,
    sexo,
    raca,
    confirmado,
    cod_ibge,
    diabetes,
    obesidade,
    imunossupressao,
    dcardiaca,
    puerpera,
    gestante,
    drc,
    assintomatico,
    tipo_teste,
    sistema
  )

#### Recode values and names to posterior bind with SG / Select only necessary columns
srag_long <- srag_long %>% ungroup() %>%   mutate(
  assintomatico = 0L,
  srag_cs_raca = as.integer(srag_cs_raca),
  srag_cs_raca = if_else(srag_cs_raca == 9L, NA_integer_, srag_cs_raca),
  tipo_teste=case_when(
    !is.na(srag_pcr_sars2)~"RT-PCR",
    !is.na(srag_an_sars2)~"Antígeno",
    srag_pcr_resul %in% c("1","2","3","5") ~"RT-PCR",
    srag_res_an %in% c("1","2","3","5") ~"Antígeno",
    (srag_pcr_resul=="4" | srag_pcr_resul=="9") |(srag_an_sars2=="4" | srag_an_sars2=="9")~"Sem teste"
  )
) %>% mutate(srag_dt_coleta=case_when(
  !is.na(srag_dt_coleta)~srag_dt_coleta,
  is.na(srag_dt_coleta) & tipo_teste %in% c("Antígeno","RT-PCR") & !is.na(srag_dt_pcr)~srag_dt_pcr,
  is.na(srag_dt_coleta) & tipo_teste %in% c("Antígeno","RT-PCR") & !is.na(srag_dt_res_an)~srag_dt_res_an 
)) %>% 
  mutate(srag_confirmado=case_when((srag_pcr_sars2=="1" | srag_an_sars2=="1")~"Confirmado",
                                   datediff(srag_dt_coleta,srag_dt_sin_pri)<11 &
                                     (is.na(srag_pcr_sars2) & is.na(srag_an_sars2)) &
                                     (srag_pcr_resul=="2" | srag_an_sars2=="2") ~"Negativo")) %>% 
  filter(!is.na(srag_confirmado)) %>% 
  mutate(across(c(srag_diabetes, srag_imunodepre, srag_renal, srag_obesidade, srag_cardiopati, srag_puerpera), ~ if_else(.x==1, 1L, 0L)),
         srag_cs_gestant = case_when(
           srag_cs_gestant %in% c(1, 2, 3, 4) ~ 1L,
           srag_cs_gestant %in% c(5, 6) ~ 0L,
           TRUE ~ 0L
         )
  ) %>%
  rename(
    diabetes = srag_diabetes,
    obesidade = srag_obesidade,
    imunossupressao = srag_imunodepre,
    dcardiaca = srag_cardiopati,
    puerpera = srag_puerpera,
    drc = srag_renal,
    gestante = srag_cs_gestant,
    idade = srag_idade_srag,
    raca = srag_cs_raca,
    sexo = srag_cs_sexo,
    dt_coleta = srag_dt_coleta,
    dt_inicio_sintomas = srag_dt_sin_pri,
    dt_notificacao = srag_dt_notif,
    confirmado = srag_confirmado,
    cod_ibge = srag_co_mun_res
  ) %>%
  mutate(sistema="srag") %>% # adicionar variável indicando o tipo de notificação
  select(
    id_vigvac,
    dt_notificacao,
    dt_inicio_sintomas,
    dt_coleta,
    idade,
    sexo,
    raca,
    confirmado,
    cod_ibge,
    diabetes,
    obesidade,
    imunossupressao,
    dcardiaca,
    puerpera,
    gestante,
    drc,
    assintomatico,
    srag_hospital,
    srag_dt_interna,
    srag_evolucao,
    srag_dt_evoluca,
    tipo_teste,
    sistema)

banco_tnd <- sdf_bind_rows(srag_long,sg_long)
#spark_write_parquet(banco_tnd,"/dados/Analysis/Thiago/Scripts/CoronaVac-Waning/db/tnd_raw")

##
banco_tnd_dt <- banco_tnd %>% select(id_vigvac:sexo,confirmado,cod_ibge,assintomatico:sistema,gestante,puerpera)
banco_tnd_clinica <- banco_tnd %>% select(id_vigvac,diabetes:dcardiaca,drc) %>% 
  group_by(id_vigvac) %>% 
  summarise(across(diabetes:drc, sum)) %>% 
  mutate(across(diabetes:drc,~if_else(.x>0,1L,0L)))


banco_raca <- banco_tnd %>%
  select(id_vigvac, dt_notificacao, raca) %>%
  dbplyr::window_order(dt_notificacao) %>%
  filter(!is.na(raca)) %>% 
  group_by(id_vigvac) %>% 
  mutate(valid=row_number()) %>%
  filter(valid==1) %>% 
  select(-c(valid,dt_notificacao)) %>% 
  ungroup()

banco_tnd_pos <- banco_tnd_dt %>% 
  left_join(banco_tnd_clinica,by="id_vigvac") %>% 
  left_join(banco_raca,by="id_vigvac")

banco_tnd_pos <- banco_tnd_pos %>% mutate(across(gestante:drc, ~ if_else(is.na(.x), 0L, .x)))
spark_write_parquet(banco_tnd_pos,"/dados/Analysis/Thiago/DB_processed/tnd_processed111121")


