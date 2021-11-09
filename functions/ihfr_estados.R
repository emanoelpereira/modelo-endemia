if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)}
if(!require(lubridate)){install.packages("lubridate"); library(lubridate)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
source("../nowcasting/fct/get.last.date.R")
source("../nowcasting/fct/read.sivep.R")

#' Função que calcula o IHFR para cada estado
#'
#' @param estado # Qual estado calcular a IHFR default Brasil todo 'estado'.
#' @param by_week # Se deve haver agregação por semana também 'by_week'.
#' @param more_age_bins # Mais faixas etárias, default de 20 em 20 anos 'more_age_bins'.
#' @param clean_environment # Op??o de deletar dados_brutos ap?s filtrar os dados, Default TRUE
#' @param overwrite # Carrega novamente os dados da SIVEP mesmo que j? exista no environment. Default TRUE
#'
#' @return # Data.frame IHFR por faixas etárias por estado por tipo de SRAG.
#' @importFrom magrittr %>%
#' 
ihfr_estados <- function(estado = NULL, 
                         by_week = FALSE, 
                         more_age_bins = FALSE,
                         clean_environment = TRUE,
                         overwrite = TRUE) {

#1.  Leitura dos dados ####
data.dir <- "../dados/SIVEP-Gripe/"
data.sivep <- get.last.date(data.dir)

if(!exists("dados_brutos",mode="list") | overwrite){
  dados_brutos <- read.sivep(dir = data.dir, escala = "pais", data = data.sivep)
}

if  (!is.null(estado)){
  dados_brutos <- dados_brutos %>% 
    filter(sg_uf == estado) %>%
    select(sg_uf, dt_sin_pri, nu_idade_n, hospital, pcr_sars2, classi_fin,
           evolucao)
}

######colocando em classes etárias########
dados_brutos$nu_idade_n <- as.numeric(dados_brutos$nu_idade_n)

if (more_age_bins == TRUE){
  dados <- dados_brutos  %>%
      mutate(age_clas = case_when(nu_idade_n >= 0 & nu_idade_n <= 4   ~ 0,
                                  nu_idade_n >= 5 & nu_idade_n <= 9   ~ 5,
                                  nu_idade_n >= 10 & nu_idade_n <= 14 ~ 10,
                                  nu_idade_n >= 15 & nu_idade_n <= 19 ~ 15,
                                  nu_idade_n >= 20 & nu_idade_n <= 24 ~ 20,
                                  nu_idade_n >= 25 & nu_idade_n <= 29 ~ 25,
                                  nu_idade_n >= 30 & nu_idade_n <= 34 ~ 30,
                                  nu_idade_n >= 35 & nu_idade_n <= 39 ~ 35,
                                  nu_idade_n >= 40 & nu_idade_n <= 44 ~ 40,
                                  nu_idade_n >= 45 & nu_idade_n <= 49 ~ 45,
                                  nu_idade_n >= 50 & nu_idade_n <= 54 ~ 50,
                                  nu_idade_n >= 55 & nu_idade_n <= 59 ~ 55,
                                  nu_idade_n >= 60 & nu_idade_n <= 64 ~ 60,
                                  nu_idade_n >= 65 & nu_idade_n <= 69 ~ 65,
                                  nu_idade_n >= 70 & nu_idade_n <= 74 ~ 70,
                                  nu_idade_n >= 75 & nu_idade_n <= 79 ~ 75,
                                  nu_idade_n >= 80 & nu_idade_n <= 84 ~ 80,
                                  nu_idade_n >= 85 & nu_idade_n <= 89 ~ 85,
                                  nu_idade_n >= 90 ~ 90))
}
else {
  dados <- dados_brutos  %>%
    mutate(age_clas = case_when(nu_idade_n >= 0 & nu_idade_n <= 19 ~ "age_0_19",
                                nu_idade_n >= 20 & nu_idade_n <= 59 ~ "age_20_59",
                                nu_idade_n >= 60 ~ "age_60+"))
}

# Deixando o environment mais leve
if(clean_environment) {
  rm(dados_brutos)
  gc()
}

###############COVID##########################
covid <- dados %>%
  filter(hospital == 1) %>%
  filter(pcr_sars2 == 1 | classi_fin == 5)  %>%
  filter(evolucao == 1 | evolucao == 2) %>% 
  filter(!is.na(age_clas)) %>%
  select(dt_sin_pri, age_clas, sg_uf, evolucao)

# Classificando semana epidemiologica por estado
## Semana epidemiologica brasileira
covid$week <- epiweek(covid$dt_sin_pri) ####semana epidemiologica começando no domingo

###################SRAG#################################
srag <- dados %>%
  filter(hospital == 1) %>%
  filter(evolucao == 1 | evolucao == 2) %>% 
  filter(!is.na(age_clas)) %>%
  select(dt_sin_pri, age_clas,sg_uf, evolucao)

# Classificando semana epidemiologica por estado
## Semana epidemiologica brasileira
srag$week <- epiweek(srag$dt_sin_pri) ####semana epidemiol?gica come?ando no domingo

if (by_week == TRUE){
  ##COVID##
  covid_ihfr  <-
    covid %>%
    dplyr::group_by(week, age_clas, sg_uf) %>%
    dplyr::summarise(sobre = sum(evolucao == 1), obitos = sum(evolucao == 2))
  
  covid_ihfr$ihfr <- covid_ihfr$obitos/(covid_ihfr$obitos + covid_ihfr$sobre)
  
  ##SRAG###
  srag_ihfr  <-
    srag %>%
    dplyr::group_by(week, age_clas, sg_uf) %>%
    dplyr::summarise(sobre = sum(evolucao == 1), obitos = sum(evolucao == 2))
  
  srag_ihfr$ihfr <- srag_ihfr$obitos/(srag_ihfr$obitos + srag_ihfr$sobre)
  
  ihfr_states <- merge(srag_ihfr, covid_ihfr, 
                       by = c("week", "age_clas", "sg_uf"),
                       suffixes = c(".srag", ".covid")) %>% 
    select(week, sg_uf, age_clas, ihfr.srag, ihfr.covid) %>% 
    as.data.frame()  %>% arrange(week)
  return(ihfr_states)
}

##COVID##
covid_ihfr  <-
  covid %>%
  dplyr::group_by(age_clas, sg_uf) %>%
  dplyr::summarise(sobre = sum(evolucao == 1), obitos = sum(evolucao == 2))

covid_ihfr$ihfr <- covid_ihfr$obitos/(covid_ihfr$obitos + covid_ihfr$sobre)

##SRAG###
srag_ihfr  <-
  srag %>%
  dplyr::group_by(age_clas, sg_uf) %>%
  dplyr::summarise(sobre = sum(evolucao == 1), obitos = sum(evolucao == 2))

srag_ihfr$ihfr <- srag_ihfr$obitos/(srag_ihfr$obitos + srag_ihfr$sobre)

ihfr_states <- merge(srag_ihfr, covid_ihfr, 
                     by = c("age_clas", "sg_uf"),
                     suffixes = c(".srag", ".covid")) %>% 
  dplyr::select(sg_uf, age_clas, ihfr.srag, ihfr.covid) %>% 
  as.data.frame()

return(ihfr_states)

}
# ####################summary dos casos UTI ########################
# 
# ##Casos##
# ##COVID##
# tabela_uti  <-
#   covid %>%
#   filter(uti == 1 | uti == 2) %>%
#   group_by(week, age_clas, sg_uf) %>%
#   summarise(N = n(),
#             leito = sum(uti == 2), uti = sum(uti == 1))
# 
# tabela_uti$uti_obs <- tabela_uti$uti / (tabela_uti$leito + tabela_uti$uti)
# 
# ##SRAG###
# tabela2_uti  <-
#   srag %>%
#   filter(uti == 1 | uti == 2) %>%
#   group_by(week, age_clas, sg_uf) %>%
#   summarise(N = n(),
#             leito = sum(uti == 2), uti = sum(uti == 1))
# 
# tabela2_uti$uti_obs <- tabela2_uti$uti/(tabela2_uti$leito + tabela2_uti$uti)
# 
# write.csv(tabela_uti, file = paste0("output/dados/", "summary_covid_IFHR_uti","_", last.date,".csv"),
#           row.names = FALSE)
# write.csv(tabela2_uti, file = paste0("output/dados/", "summary_srag_IFHR_uti","_", last.date,".csv"),
#           row.names = FALSE)
# 
# ##Death##
# ##COVID##
# tabela_uti_death  <-
#   covid %>%
#   filter(uti == 1 | uti == 2) %>%
#   filter(evolucao == 2) %>%
#   group_by(week, age_clas, sg_uf) %>%
#   summarise(N = n(),
#             leito = sum(uti == 2), uti = sum(uti == 1))
# 
# tabela_uti_death$uti_obs <- tabela_uti_death$uti / (tabela_uti_death$leito + tabela_uti_death$uti)
# 
# ##SRAG###
# tabela2_uti_death  <-
#   srag %>%
#   filter(uti == 1 | uti == 2) %>%
#   filter(evolucao == 2) %>%
#   group_by(week, age_clas, sg_uf) %>%
#   summarise(N = n(),
#             leito = sum(uti == 2), uti = sum(uti == 1))
# 
# tabela2_uti_death$uti_obs <- tabela2_uti_death$uti / (tabela2_uti_death$leito + tabela2_uti_death$uti)
# 
# 
# write.csv(tabela_uti_death, file = paste0("output/dados/", "summary_covid_IFHR_uti_obs","_", last.date,".csv"),
#           row.names = FALSE)
# write.csv(tabela2_uti_death, file = paste0("output/dados/", "summary_srag_IFHR_uti_obs","_", last.date,".csv"),
#           row.names = FALSE)
# 
# ####HOSPITALIZADOS
# 
# hosp_week <- srag %>%
#   group_by(week, sg_uf) %>%
#   summarise(hosp = n())
# 
# hosp_week <- hosp_week %>%
#   filter(week < 36 & week >= 10)
# 
# write.csv(hosp_week, file = paste0("output/dados/", "hospitalizados_srag_week","_", last.date,".csv"),
#           row.names = FALSE)
