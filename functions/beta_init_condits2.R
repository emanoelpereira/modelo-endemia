#' beta_init_condits2.
#' 
#' A function that calculates the initial conditions and the average contact rate and success between unvaccinated susceptibles and infectious individuals
#' given the age distributed new hospitalizations by week. The method is based on the linearization of the system.
#'
#' @param r # Fitted exponetial growth rate in 1/days 'r'.
#' @param NEW.HOSP.AGE # New hospitalizations per week by age, a vector with 3 elements: HJ, HA, HO 'NEW.HOSP.AGE'.
#' @param PREVALENCE # Fraction of the population that already had the disease by age, a vector with 3 elements: PREVJ, PREVA, PREVO 'PREVALENCE'.
#' @param POP.DISTR # Population distributed by age, a vector with 3 elements: POPJ, POPA, POPO 'POP.DISTR'.
#' @param CONTACT.M # A contact matrix, must give as matrix 'CONTACT.M'.
#' @param EXPOSURE.PERIOD.DAYS # Average time between being infected and developing symptoms 'EXPOSURE.PERIOD.DAYS'.
#' @param SICKNESS.PERIOD.DAYS # Average time between being infectious and recovering for asymptomatic and mild 'SICKNESS.PERIOD.DAYS'.
#' @param SEVERE.PERIOD.DAYS # Average time between being infectious and recovering/dying for severe cases 'SEVERE.PERIOD.DAYS'.
#' @param CONT.REDUC.FRAC # Reduction on the expose of symptomatic (due to symptoms/quarantining) 'CONT.REDUC.FRAC'.
#' @param SEVERE.CONT.REDUC.FRAC # Reduction on the expose of severe cases (due to hospitalization) 'SEVERE.CONT.REDUC.FRAC'. 
#'                                 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
#' @param REL.INFEC.PRESYMP # relative infectiousness of pre-symptomatic individuals 'REL.INFEC.PRESYMP'.
#' @param ASYMPTOMATIC.FRAC # Fraction of asymptomatic cases in total cases 'ASYMPTOMATIC.FRAC'.
#' @param SEVERITY.FRAC # Fraction of severe cases/hospitalizations in symptomatic cases (IHR) 'SEVERITY.FRAC'.
#' @param DEATH.FRAC # Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) 'DEATH.FRAC'.
#' @param V2.FRAC # Fraction of the infected people due to the second strain 'V2.FRAC'.

#'
#' @return A list where the first entry is BETA.RATE the average contact rate and success between unvaccinated susceptibles and infectious individuals
#'         and the second entry is the vector with the initial conditions.
#' @export
#' @import matlib, rARPACK
#'
#' @examples

library(rARPACK)
library(matlib)

init_condits <- function(r, NEW.HOSP.AGE, PREVALENCE, POP.DISTR, CONTACT.M, EXPOSURE.PERIOD.DAYS,
                         SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                         SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                         SEVERITY.FRAC, DEATH.FRAC, V2.FRAC = 0.0001){

  EXPOSURE.PERIOD.WEEKS <- EXPOSURE.PERIOD.DAYS / 7.0
  SICKNESS.PERIOD.WEEKS <- SICKNESS.PERIOD.DAYS / 7.0
  SEVERE.PERIOD.WEEKS <- SEVERE.PERIOD.DAYS / 7.0
  r <- unname(r) # measured in 1/day
  
  POP.C1 <- c(0.0, 0.0, 0.0)
  POP.R1 <- POP.DISTR * PREVALENCE
  POP.D1 <- c(0.0, 0.0, 0.0)

  POP.C2 <- c(0.0, 0.0, 0.0)
  POP.R2 <- c(0.0, 0.0, 0.0)
  POP.D2 <- c(0.0, 0.0, 0.0)

  POP.S <- POP.DISTR - POP.R1
  
  BETA.RATE <- bissec(function(BETA.RATE) beta_relation_r(BETA.RATE, r, POP.S, sum(POP.DISTR), CONTACT.M, EXPOSURE.PERIOD.DAYS,
                                                           SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                                                           SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                                                           SEVERITY.FRAC, DEATH.FRAC), 0.0, 0.1)
  
  POP.AUTO <- r_related_eigen(BETA.RATE, POP.S, sum(POP.DISTR), CONTACT.M, EXPOSURE.PERIOD.DAYS,
                               SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                               SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                               SEVERITY.FRAC, DEATH.FRAC)

  # TODO: faz diferença normalizar por AUTO[1] ou AUTO[3], pois os valores podem ser diferentes! Como fazer para garantir os valores de E sejam os abaixo?
  #POP.AUTO <- Re(POP.AUTO / POP.AUTO[1])
  # normaliza cada idade pela entrada daquela idade?!
  POP.AUTO <- Re(POP.AUTO / POP.AUTO[1:3])

  POP.E1 <- (1 - V2.FRAC) * EXPOSURE.PERIOD.WEEKS * NEW.HOSP.AGE / SEVERITY.FRAC
  POP.E2 <- V2.FRAC * EXPOSURE.PERIOD.WEEKS * NEW.HOSP.AGE / SEVERITY.FRAC

  POP.INF1 <- POP.E1 * POP.AUTO
  POP.INF2 <- POP.E2 * POP.AUTO
  # Colocar sum gera problemas de população negativa em idosos
  POP.S <- POP.DISTR - (POP.INF1[1:3] + POP.INF1[4:6] + POP.INF1[7:9] + POP.INF1[10:12] + POP.INF2[1:3] + POP.INF2[4:6] + POP.INF2[7:9] + POP.INF2[10:12] + POP.D1 + POP.R1 + POP.R2 + POP.D2)

  POP0 <- c(POP.S, POP.INF1, POP.C1, POP.R1, POP.D1, 
                   POP.INF2, POP.C2, POP.R2, POP.D2)
  
  return(list(BETA.RATE = BETA.RATE, POP0 = POP0))
}

beta_relation_r <- function(BETA.RATE, rfit, S, P, CONTACT.M, EXPOSURE.PERIOD.DAYS,
                            SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                            SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                            SEVERITY.FRAC, DEATH.FRAC){
  #following the E,A,I,H order
  sml = t(t(CONTACT.M)%*%diag(x=S))
  f = cbind(REL.INFEC.PRESYMP*sml,
                 sml,
                 (1-CONT.REDUC.FRAC)*sml,
                 (1-SEVERE.CONT.REDUC.FRAC)*sml)
  f_full = BETA.RATE/P*rbind(f,matrix(0, nrow = 9, ncol =12))
  unit <- c(1,1,1)
  v1 = cbind(diag(1/EXPOSURE.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=9))
  v2 = cbind(diag(-ASYMPTOMATIC.FRAC*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), diag(1/SICKNESS.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=6))
  v3 = cbind(diag(-(1 - ASYMPTOMATIC.FRAC)*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3),
            diag(1/SICKNESS.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3))
  v4 = cbind(diag(-SEVERITY.FRAC/EXPOSURE.PERIOD.DAYS*unit),
           matrix(0,nrow=3,ncol=6),
           diag((1/SEVERE.PERIOD.DAYS)*unit))
  V = rbind(v1,v2,v3,v4)

  values <- eigs(f_full - V, 1, which = 'LR')$values
  return(Re(values) - rfit)
}

r_related_eigen <- function(BETA.RATE, S, P, CONTACT.M, EXPOSURE.PERIOD.DAYS,
                            SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                            SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                            SEVERITY.FRAC, DEATH.FRAC){
  #following the E,A,I,H order
  sml = t(t(CONTACT.M)%*%diag(x=S))
  f = cbind(REL.INFEC.PRESYMP*sml,
                 sml,
                 (1-CONT.REDUC.FRAC)*sml,
                 (1-SEVERE.CONT.REDUC.FRAC)*sml)
  f_full = BETA.RATE/P*rbind(f,matrix(0, nrow = 9, ncol =12))
  unit <- c(1,1,1)
  v1 = cbind(diag(1/EXPOSURE.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=9))
  v2 = cbind(diag(-ASYMPTOMATIC.FRAC*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), diag(1/SICKNESS.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=6))
  v3 = cbind(diag(-(1 - ASYMPTOMATIC.FRAC)*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3),
            diag(1/SICKNESS.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3))
  v4 = cbind(diag(-SEVERITY.FRAC/EXPOSURE.PERIOD.DAYS*unit),
           matrix(0,nrow=3,ncol=6),
           diag((1/SEVERE.PERIOD.DAYS)*unit))
  V = rbind(v1,v2,v3,v4)

  vectors <- eigs(f_full - V, 1, which = 'LR')$vectors
  return(Re(vectors))
}

bissec <- function(f, a, b, n = 1000, tol = 0.000000001){
  i = 0
  while(i<n){
    if(f(a)==0.0){
      result <- a
      break
    }
    if(f(b)==0.0){
      result <- b
      break
    }
    c <- (a + b)/2.0
    if(f(c)==0|(b-c)<tol){
      result <- c
      break
    }
    i <- i + 1
    if(sign(f(c))==sign(f(a))){return(bissec(f, c, b))}
    else{return(bissec(f, a, c))}
  }
  if (!result){print('Root finding failed')}
  else{return(result)}
}
