# logit and inverse logit functions
logit <- function(x) log(x/(1-x))
inv.logit <- function(x) exp(x)/(1+exp(x))

