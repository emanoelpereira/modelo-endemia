library(dplyr)

if(!exists("aggregate_age_matrix",mode="function")) source('./functions/aggregate_age_matrix.R')
if(!exists("carregar_dados_estado",mode="function")) source('./functions/dados_estados.R')
if(!exists("init_condits",mode="function")) source('./functions/beta_init_condits2.R')
if(!exists("fitP.exp",mode="function"))  source('./functions/fitP.exp.R')

######################################
##### TWO-STRAIN SEIRHD MODEL ########
######################################
# This file contains the epidemiological
# parameters for COVID-19 modelling in
# Brazilian population.

#####################################
######## General Notation ###########
#####################################
# About the units: .DAYS in days
#                  .RATE in /day
#                  .FRAC in [0, 1]
#                  .NUM  integer
# About the ages:
# 3 age classes: juvenile (J) <20yrs
#                adults (A) 20-60yrs
#                elderly(I)   >60yrs
# Parameters which depend on age are
# represented as vectors in order: J A I


# Parameters of the initial population

# Population Age Distribution
all.age.distr <- read.csv('DATA/DistrEtaria2020.csv')
if (! exists("Estado")) {
    Estado  <- "SP"
}
age.distr <- all.age.distr[-1, Estado]
POP.TOTAL.NUM = sum(age.distr) # Total Population Size
POP.DISTR <- c(sum(age.distr[1:4]), sum(age.distr[5:12]), sum(age.distr[13:length(age.distr)]))

# Contact matrix    ( JJ, JA, JI, AJ, AA, AI, IJ, IA, II)
#CONTACT.M = matrix(c(1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3), nrow=3)
load('DATA/prem2020_bra.Rdata')

full.mat.cont <- new_all
# alternate matrix considering different internventions
# full.mat.cont <- new_home + 0.7*new_others + 0.7*new_work + 0.2*new_school

# aggregated Contact Matrix
aggregate_indices <- list(c(1:4), c(5:12), c(13:16))
CONTACT.M <- aggregate_age_matrix(mat.cont = full.mat.cont,
                                  aggregate_indices = aggregate_indices,
                                  age_structure = age.distr)

# the names with a preceding \ indicate the parameter name in the model equations
# the names in () represent the respective parameter name in CoMo BR model
# references are in the CoMo BR model parameter table

EXPOSURE.PERIOD.DAYS   = 5.8  # (gamma) Average time between being infected and developing symptoms \gamma
SICKNESS.PERIOD.DAYS   = 9    # (nui) Average time between being infectious and recovering for asymptomatic and mild \nu_i
SEVERE.PERIOD.DAYS     = 8.4  # (nus) Average time between being infectious and recovering/dying for severe cases \nu_s

CONT.REDUC.FRAC        = 0.1  # Reduction on the expose of symptomatic (due to symptoms/quarantining). \tau 0 means the same level of exposure of asymptomatic and 1 means no expose whatsoever. 
SEVERE.CONT.REDUC.FRAC = 0.9  # Reduction on the expose of severe cases (due to hospitalization). \eta 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
REL.INFEC.PRESYMP      = 1    # (rho) relative infectiousness of pre-symptomatic individuals. \rho

#### Classification of cases related to the severity of cases of unvaccinated individuals
ASYMPTOMATIC.FRAC = c(0.67, 0.44, 0.31) # Fraction of asymptomatic cases in total cases (pclin) \alpha

# Fraction of severe cases/hospitalizations in symptomatic cases (IHR) \sigma
ihr <- read.csv('DATA/ihr.csv')
weighted_ihr <- ihr[,2]/100 * age.distr
SEVERITY.FRAC = c(sum(weighted_ihr[1:4]) / sum(age.distr[1:4]),
                  sum(weighted_ihr[5:12]) / sum(age.distr[5:12]),
                  sum(weighted_ihr[13:19]) / sum(age.distr[13:19]))

# Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) \mu
ihfr_estados <- read.csv('DATA/ihfr.csv')
DEATH.FRAC <- ihfr_estados %>%
    filter(sg_uf == Estado) %>%
    arrange(age_clas) %>%
    .$ihfr.covid

# probability of reporting hospitalized case
# (could be, for instance, ratio of Covid / SRAG reported)
REPORT.PROB = 1.

### NEW VARIANT parameters
REL.BETA.RATE2 = 1.3
REL.SEVERITY.FRAC.V2 = 1.0
PROB.REINFEC.V2 = 0.
EXPOSURE.PERIOD.DAYS.V2   = 5.8  # (gamma) Average time between being infected and developing symptoms \gamma
SICKNESS.PERIOD.DAYS.V2   = 9    # (nui) Average time between being infectious and recovering for asymptomatic and mild \nu_i
SEVERE.PERIOD.DAYS.V2     = 8.4  # (nus) Average time between being infectious and recovering/dying for severe cases \nu_s

CONT.REDUC.FRAC.V2        = 0.1  # Reduction on the expose of symptomatic (due to symptoms/quarantining). \tau 0 means the same level of exposure of asymptomatic and 1 means no expose whatsoever. 
SEVERE.CONT.REDUC.FRAC.V2 = 0.9  # Reduction on the expose of severe cases (due to hospitalization). \eta 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
REL.INFEC.PRESYMP.V2      = 1    # (rho) relative infectiousness of pre-symptomatic individuals. \rho
ASYMPTOMATIC.FRAC.V2 = ASYMPTOMATIC.FRAC
DEATH.FRAC.V2 = DEATH.FRAC
SEVERITY.FRAC.V2 = REL.SEVERITY.FRAC.V2*SEVERITY.FRAC

### Initial conditions

## to use updated data from SIVEP: modify and run functions/update_data_from_sivep.R
## exports cases_dist_age, recent_cases_dist_age, and ihfr

# prevalence
cases_dist_age <- read.csv('DATA/cases_dist_age_states.csv')
cases_dist_age <- cases_dist_age %>%
    filter(sg_uf == Estado) %>%
    arrange(age_clas) %>%
    .$N.covid
Recovered <- cases_dist_age / (ihr[,2] / 100)
POP.R1 <- c(sum(Recovered[1:4]),
            sum(Recovered[5:12]),
            sum(Recovered[13:19]))

# age distribution of recent hospitalized cases
recent_cases_dist_age <- read.csv('DATA/recent_cases_dist_age_states.csv')
recent_cases_dist_age <- recent_cases_dist_age %>%
    filter(sg_uf == Estado) %>%
    arrange(age_clas) %>%
    .$N.covid

# Load srag and covid data from specified state
state.zoo.list <- carregar_dados_estado(estado = Estado, data.base = "2021_02_15")
#srag.zoo <- diff(state.zoo.list$now.srag.zoo)
covid.zoo <- diff(state.zoo.list$now.covid.zoo)
data.base <- state.zoo.list$data.base
# exponential fit
r <- fitP.exp(tail(covid.zoo, n = 3), only.coef = FALSE)$coefficients[2]

new.hosp <- coredata(covid.zoo[max(time(covid.zoo))]) * recent_cases_dist_age / sum(recent_cases_dist_age)

init.conds <- init_condits(r, new.hosp, PREVALENCE = POP.R1/POP.DISTR, POP.DISTR,
             CONTACT.M, EXPOSURE.PERIOD.DAYS, SICKNESS.PERIOD.DAYS,
             SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC, SEVERE.CONT.REDUC.FRAC,
             REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC, SEVERITY.FRAC, DEATH.FRAC,
             V2.FRAC = 0)

POP0 <- init.conds$POP0
# Average contact rate and success between unvaccinated susceptibles and
# infectious individuals \beta
# calculated from hospitalized cases in the last 3 weeks
BETA.RATE1 <- init.conds$BETA.RATE
names(POP0) <- paste0(rep(c("POP.S", "POP.E1", "POP.A1", "POP.I1", "POP.H1",
                             "POP.C1", "POP.R1", "POP.D1",
                             "POP.E2", "POP.A2", "POP.I2", "POP.H2",
                             "POP.C2", "POP.R2", "POP.D2"), each=3), rep(1:3, 15))


## Removed variables used in this code but not needed in the future
## in order to keep the workspace clear
rm(ihfr_estados, weighted_ihr, all.age.distr, age.distr, new_all, new_home,
   new_others, new_school, new_work, u_all, u_home, u_work, u_others, u_school,
   ihr, full.mat.cont, aggregate_indices, r_all, r_home, r_others,
   r_school, r_work, cases_dist_age, recent_cases_dist_age, Recovered, POP.R1,
   new.hosp)
