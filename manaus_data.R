################################################################################
##  location-specific data
################################################################################
source('functions/aggregate_age_matrix.R')
source('functions/dados_municipios.R')
source('functions/logit.R')

# dados demográficos
# fonte: https://demografiaufrn.net/laboratorios/lepp/
age.distr <- read.csv('DATA/DistrEtaria2020_Manaus.csv')
POP.TOTAL.NUM = sum(age.distr) # Total Population Size
POP.DISTR <- c(sum(age.distr[1:4]), sum(age.distr[5:12]), sum(age.distr[13:length(age.distr)]))

# Contact matrix    ( JJ, JA, JI, AJ, AA, AI, IJ, IA, II)
load('DATA/prem2020_bra.Rdata')

full.mat.cont <- new_all
# alternate matrix considering different internventions
# full.mat.cont <- new_home + 0.7*new_others + 0.7*new_work + 0.2*new_school

# aggregated Contact Matrix
aggregate_indices <- list(c(1:4), c(5:12), c(13:16))
CONTACT.M <- aggregate_age_matrix(mat.cont = full.mat.cont,
                                  aggregate_indices = aggregate_indices,
                                  age_structure = age.distr)

# Fraction of severe cases/hospitalizations in symptomatic cases (IHR) \sigma
ihr <- read.csv('DATA/ihr.csv')
weighted_ihr <- ihr[,2]/100 * age.distr
SEVERITY.FRAC = c(sum(weighted_ihr[1:4]) / sum(age.distr[1:4]),
                  sum(weighted_ihr[5:12]) / sum(age.distr[5:12]),
                  sum(weighted_ihr[13:19]) / sum(age.distr[13:19]))

# Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) \mu
# usando do AM mesmo, deve ter pouca diferença
ihfr_estados <- read.csv('DATA/ihfr.csv')
DEATH.FRAC <- ihfr_estados %>%
    filter(sg_uf == "AM") %>%
    arrange(age_clas) %>%
    .$ihfr.covid

# age distribution of recent hospitalized cases
# usando do AM mesmo, deve ter pouca variação
recent_cases_dist_age <- read.csv('DATA/recent_cases_dist_age_states.csv')
recent_cases_dist_age <- recent_cases_dist_age %>%
    filter(sg_uf == "AM") %>%
    arrange(age_clas) %>%
    .$N.covid

# Load srag and covid data from specified state
state.zoo.list <- carregar_dados_cidades(estado = "AM", cidade = "Manaus")
#srag.zoo <- diff(state.zoo.list$now.srag.zoo)
covid.zoo <- diff(state.zoo.list$now.covid.zoo)
data.base <- state.zoo.list$data.base

new.hosp <- as.numeric(covid.zoo[as.Date("2020-11-01")]) *
    recent_cases_dist_age / sum(recent_cases_dist_age)

