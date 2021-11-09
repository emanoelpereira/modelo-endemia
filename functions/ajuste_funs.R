library(doParallel)
library(dplyr)
library(lubridate)
library(minpack.lm)
library(zoo)
library(rARPACK)
library(pse)
library(bbmle)
library(gridExtra)
source('SEIR_run.R')
source('./functions/beta_init_condits2.R')
source('functions/logit.R')
source('functions/end.of.epiweek.R')
source('functions/C_calculator.R')

model_solution <- function(params, full_solution = FALSE, ...) {
    if ("prevalence" %in% names(params) | "r" %in% names(params)) {
        if ("prevalence" %in% names(params))
            PREVALENCE <- params["prevalence"] / 100 * c(1,1,1)
        if ("r" %in% names(params))
            r <- params["r"]
        # adjust initial conditions
        init.conds <- init_condits(r, new.hosp, PREVALENCE = PREVALENCE, POP.DISTR,
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
    }
    if ("IHR.V2.PROP" %in% names(params)) {
        ihr2 <- ihr[,2]/100
        ihr2 <- 1 / (1 + 1/(params["IHR.V2.PROP"] * ihr2/(1-ihr2)))
        weighted_ihr2 <- ihr2 * age.distr
        SEVERITY.FRAC.V2 = c(sum(weighted_ihr2[1:4]) / sum(age.distr[1:4]),
                             sum(weighted_ihr2[5:12]) / sum(age.distr[5:12]),
                             sum(weighted_ihr2[13:length(age.distr)]) / sum(age.distr[13:length(age.distr)]))
    }

    # introduce new variant in initial conditions
    new.POP0 <- POP0
    new.POP0[4:15] <- POP0[4:15] * (1-params["INIT.V2.FRAC"])
    new.POP0[25:36] <- POP0[4:15] * params["INIT.V2.FRAC"]

    diffEqs <- func.factory(REL_BETA_RATE2 = params["REL.BETA.RATE2"],
                            PROB_REINFEC_V2 = params["PROB.REINFEC.V2"],
                            BETA_RATE1 = BETA.RATE1,
                            SEVERITY_FRAC_V2 = SEVERITY.FRAC.V2,
                            ...)

    SOLUTION <- ode(y = new.POP0, times = TIME.VECTOR, func = diffEqs, parms = c())
    SOLUTION <- C_calculator(SOLUTION, EXPOSURE.PERIOD.DAYS.V2, SEVERITY.FRAC.V2)
    if (full_solution)
        return(SOLUTION)
    sol <- data.frame(time = SOLUTION[-1,"time"],
                      C1 = diff(SOLUTION[,"POP.C11"] + SOLUTION[,"POP.C12"] + SOLUTION[,"POP.C13"]),
                      C2 = diff(SOLUTION[,"POP.C21"] + SOLUTION[,"POP.C22"] + SOLUTION[,"POP.C23"]))
    sol$time <- d0 + (sol$time -1)

    # aggregating data by the frequency time windows
    freq <- rep(NA, nrow(d))
    for (i in 1:nrow(d)) {
        soli <- sol[sol$time >= d[i,"t_start"] & sol$time <= d[i,"t_end"],]
        freq[i] <- sum(soli$C2) / sum(soli$C1, soli$C2)
    }

    # aggregating by epi week
    sol <- sol %>%
        mutate(week = end.of.epiweek(time, end = 0)) %>%
        group_by(week) %>%
        summarise(C1 = sum(C1), C2 = sum(C2)) %>%
        as.data.frame()
    sol.zoo <- zoo(sol[,c("C1","C2")], sol$week)
    # first point is not a full week
    sol.zoo <- sol.zoo[-1]

    d.cases <- merge.zoo(sol.zoo, covid.zoo, all = FALSE)
    return(list(cases = d.cases, freq = cbind(d, pred.freq=freq)))
}

residual_sol <- function(sol, weight.freq = 1, use.logit = TRUE) {
    cases <- sol$cases
    freq <- sol$freq$pred.freq
    cases <- (cases - mean(cases$covid.zoo)) / sd(cases$covid.zoo)
    cases <- cases$C1 + cases$C2 - cases$covid.zoo
    if (use.logit) {
        freq.logit <- logit(freq)
        if (any(! is.finite(c(freq.logit, d$freq.logit))))
            warning("Infinite logit value!")
        d.freq <- (freq.logit - d$freq.logit) / sd(d$freq.logit)
    } else {
        d.freq <- (freq - d$freq) / sd(d$freq - mean(d$freq))
    }
    return(c(as.numeric(cases), weight.freq * as.numeric(d.freq)))
}

neg_loglike_sol <- function(sol) {
    cases <- sol$cases
    freq <- sol$freq$pred.freq
    return(-(sum(dpois(cases$covid.zoo, cases$C1 + cases$C2, log = TRUE)) +
             sum(dbinom(x = d$P.1, size = d$total, prob = freq, log = TRUE))))
}

all.params.exploration <- function(lower, upper, conv.guess, LHS = TRUE, n = 100, ...) {
    if (LHS) {
        # use LHS to sample the parameter space
        q <- rep('qunif', length(lower))
        q.args <- list()
        for (p in names(lower))
            q.args[[p]] <- list(min=unname(lower[p]), max=unname(upper[p]))
        uncoupledLHS <- LHS(model=NULL, factors=names(lower), N=n,
                            q=q, q.arg=q.args, method='random')
        starters <- uncoupledLHS$data
    } else {
        N <- ceiling(n**(1/length(lower)))
        arg.ranges <- list()
        for (i in 1:length(lower))
            arg.ranges[[i]] <- seq(lower[i], upper[i], length.out = N)
        starters <- do.call(expand.grid, arg.ranges)
        colnames(starters) <- names(lower)
    }

    if (! missing(conv.guess)) {
        starters <- conv.guess(starters)
    }

    ExtraArgs <- list(...)
    if ("IHR.V2.PROP" %in% names(ExtraArgs)) {
        starters$IHR.V2.PROP  <- ExtraArgs[["IHR.V2.PROP"]]
        ExtraArgs["IHR.V2.PROP"] <- NULL
    }
    if ("prevalence" %in% names(ExtraArgs)) {
        starters$prevalence <- ExtraArgs[["prevalence"]]
        ExtraArgs["prevalence"] <- NULL
    }

    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(starters), .combine = rbind) %dopar% {
        sol <- do.call(model_solution, c(list(params=unlist(starters[i,])),
                                         ExtraArgs))
        #residual <- residual_sol(sol)
        neg.loglike <- neg_loglike_sol(sol)

        res <- c(unlist(starters[i,]), #ss = sum(residual**2),
                 neg.loglike = neg.loglike)
    }

    results <- as.data.frame(results)
}

helper.nls.fun <- function(...) {
    return(residual(model_solution(...)))
}

fit.nls.guesses <- function(guesses, lower, upper, control, ...) {
    if (missing(control)) {
        control <- nls.lm.control(ptol=1.490116e-12, ftol=1.490116e-12, maxiter=200,
                          maxfev=4*100*(length(lower) + 1))
    }
    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(guesses), .combine = rbind) %dopar% {
        parms <- unlist(guesses[i,])
        fitval <- nls.lm(par = parms, fn = helper.nls.fun, upper = upper,
                         lower = lower, control = control, ...)
        residuo <- sum(residuals(fitval)*residuals(fitval))
        neg.likelihood <- neg_loglike_sol(model_solution(coef(fitavl)))
        res <- c(coef(fitval), residuo = residuo, neg.loglike = neg.likelihood)
    }
    return(as.data.frame(results))
}

fit.mle.guesses <- function(guesses, helper, ...) {
    registerDoParallel(cores=detectCores()-1)
    results <- foreach(i = 1:nrow(guesses), .combine = rbind) %dopar% {
        guess <- guesses[i,]
        fit <- try({mle2(helper, start = as.list(guess), ...)})
        if(class(fit) == "try-error") {
            res <- c(guess, loglike = NA)
        } else {
            res <- c(coef(fit), loglike = logLik(fit))
        }
    }
    return(as.data.frame(results))
}

fit.all.prevalences <- function(lower, upper, prevalences, n=1000,
                                prop.fit=0.01, r=0, savefile) {
    fit.prevalence <- list()
    for (prevalence in prevalences) {
        PREVALENCE <- prevalence / 100 * c(1,1,1)
        # adjust initial conditions
        init.conds <- init_condits(r, new.hosp, PREVALENCE = PREVALENCE, POP.DISTR,
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

        # some brute force
        #TODO: this no longer works
        result = all.params.exploration(lower, upper, LHS=F, n=n, use.logit=T,
                                        POP0 = POP0, BETA_RATE1 = BETA.RATE1)
        # prune using sum of squares of residuals
        pruned <- result %>% slice_min(residuo, prop=prop.fit)
        # prune using negative log-likelihood
        #pruned <- result %>% slice_min(neg.loglike, prop=prop.fit)

        # preliminary explorations show that...
        extra.guesses <- data.frame(INIT.V2.FRAC=c(0.02, 0.08),
                                    REL.BETA.RATE2=c(2, 1.4),
                                    PROB.REINFEC.V2=c(0.24, 1))
        pruned <- rbind(pruned[,1:3], extra.guesses)
        # refine top 1% best fits
        best.fits <- fit.nls.guesses(pruned, lower, upper, POP0 = POP0,
                                     BETA_RATE1 = BETA.RATE1, use.logit = T)
        best.fits <- arrange(best.fits, neg.loglike)
        # save full info
        if (! missing(savefile)) {
            write.csv(best.fits, paste0(savefile, prevalence, '.csv'),
                      row.names = FALSE)
        }
        # eliminate very similar rows
        best.fits <- best.fits[row.names(unique(round(best.fits, 4))),]
        fit.prevalence[[as.character(prevalence)]] <- best.fits
    }
    results <- plyr::ldply(fit.prevalence, .id="prevalence")
    results$prevalence <- as.numeric(as.character(results$prevalence))
    return(results)
}

## Pi: adicionei um argumento para receber a tabela com as soluções, retornada
## pela funcao acima. Se isso falta e sao dados os parametros ai é feita a
## soluçao numerica de novo. Mas parece que esta está com algum bug (UPDATE: SOLVED!)
plot.fitted.sol <- function(params, solution.table, ...) {
    if (!missing(params)) {
        sol <- model_solution(params, ...)
        cases <- sol$cases
        freq <- sol$freq$pred.freq
    } else if (!missing(solution.table)) {
        cases  <- solution.table$cases
        freq <- solution.table$freq$pred.freq
    } else {
        stop("Provide parameters or solution")
    }

    par(mfrow=c(1,2))
    # plots casos
    plot(cases$covid.zoo, type="p", col="black", ylim=c(0, max(cases)),
         xlab="date of first symptoms", ylab="new weekly hospitalizations")
    lines(cases$C1, col="#00B4B8")
    lines(cases$C2, col="#F65D57")
    lines(cases$C1 + cases$C2, col="darkgrey")
    # plot freqs
    d$mtime <- d$t_start + difftime(d$t_end, d$t_start) / 2
    if ("source" %in% colnames(d)) {
        d.fiocruz <- d[d$source=="fiocruz",]
        d.cadde <- d[d$source=="cadde",]
        plot(d.cadde$mtime, d.cadde$freq, col = "black", ylim = c(0,1),
             xlim = c(min(d$t_start), max(d$t_end)), xlab = "date",
             ylab = "frequency of P.1")
        points(d.fiocruz$mtime, d.fiocruz$freq, col="black", pch=2)
    } else {
        plot(d$mtime, d$freq, col = "black", ylim = c(0,1),
             xlim = c(min(d$t_start), max(d$t_end)), xlab = "date",
             ylab = "frequency of P.1")
    }
    points(d$mtime, freq, col="#F65D57", pch=8)
    ##PI: to record table with predicted and observaed values
    return(invisible(list(cases=cases, freq = freq, obs.freq = d)))
}

plot.all.prevalences <- function(fits) {
    for (i in 1:nrow(fits)) {
        plot.fitted.sol(unlist(fits[i,]))
        title(main=paste("initial prevalence:", fits[i,"prevalence"], "%"))
        dev.copy(png, paste0("best_fit_prev_", fits[i,"prevalence"], ".png"), width=800, height=600)
        dev.off()
    }
}

## Plots of predicted x observed
## Helper function
plot.3 <- function(sol){
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2) %>%
        ggplot(aes(x = Index)) +
        geom_point(aes(y=covid.zoo)) +
        geom_line(aes(y = C1), col = "#00B4B8", lty =2) +
        geom_line(aes(y=C2), col = "#F65D57", lty =2) +
        geom_line(aes(y=total)) +
        theme_bw() +
        ylab("Casos hospitalizados por semana") +
        xlab("Semana")
    p2 <- sol$freq %>%
        fortify() %>%
        filter(source == "fiocruz") %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("Proporção P.1") +
        xlab("Mês")
    p3 <- sol$freq %>%
        fortify() %>%
        filter(source == "cadde") %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("Proporção P.1") +
        xlab("Dia")
    grid.arrange(p1, p2, p3, ncol =3, nrow =1)
}
plot.2 <- function(sol){
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2) %>%
        ggplot(aes(x = Index)) +
        geom_point(aes(y=covid.zoo)) +
        geom_line(aes(y = C1), col = "#00B4B8", lty =2) +
        geom_line(aes(y=C2), col = "#F65D57", lty =2) +
        geom_line(aes(y=total)) +
        theme_bw() +
        ylab("Casos hospitalizados por semana") +
        xlab("Semana")
    p2 <- sol$freq %>%
        fortify() %>%
        ggplot(aes(t_start)) +
        geom_point(aes(y = freq)) +
        geom_line(aes(y = pred.freq)) +
        theme_bw() +
        ylab("Proporção P.1") +
        xlab("Semana")
    grid.arrange(p1, p2, ncol = 2, nrow = 1)
}

# plot with confidence interval
# to produce plots for the paper, run before:
# Sys.setlocale("LC_ALL","en_US") # labels in english
# and save it using:
# dev.copy(png, "plot.CI.fit3.f2.png", width=1000, height=1800, res=250, pointsize=10)
# OR in EPS (open BEFORE plotting):
# cairo_ps("plot.CI.fit3.f2.eps", width=6, height=8, pointsize=12)
# dev.off()
plot.3.CI <- function(sol){
    tmin <- min(sol$freq$t_start[1], time(sol$cases)[1])
    tmax <- max(sol$freq$t_end[length(sol$freq$t_end)], max(time(sol$cases)))
    sol$freq$t_plot <- sol$freq$t_start + difftime(sol$freq$t_end, sol$freq$t_start) / 2
    color1 <- "blue" # #00B4B8
    color2 <- "red" # #F65D57
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2,
               total.inf = C1.q2.5 + C2.q2.5,
               total.sup = C1.q97.5 + C2.q97.5) %>%
        ggplot(aes(x = Index)) +
        geom_line(aes(y = C1), col = color1, lty =2) +
        geom_line(aes(y=C2), col = color2, lty =2) +
        geom_line(aes(y=total)) +
        geom_ribbon(aes(ymin = C1.q2.5, ymax = C1.q97.5), fill = "blue", alpha = 0.3) +
        geom_ribbon(aes(ymin = C2.q2.5, ymax = C2.q97.5), fill = "red", alpha = 0.3) +
        geom_ribbon(aes(ymin = total.inf, ymax = total.sup), fill = "black", alpha = 0.3) +
        geom_point(aes(y=covid.zoo)) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("New hospitalized cases") +
        xlab("")
    p2 <- sol$freq %>%
        fortify() %>%
        filter(source == "fiocruz") %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("")
    p3 <- sol$freq %>%
        fortify() %>%
        filter(source == "cadde") %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("date")
    grid.arrange(p1, p2, p3, ncol =1, nrow =3)
}
plot.2.CI <- function(sol){
    tmin <- min(sol$freq$t_start[1], time(sol$cases)[1])
    tmax <- max(sol$freq$t_end[length(sol$freq$t_end)], max(time(sol$cases)))
    sol$freq$t_plot <- sol$freq$t_start + difftime(sol$freq$t_end, sol$freq$t_start) / 2
    color1 <- "blue" # #00B4B8
    color2 <- "red" # #F65D57
    p1 <- sol$cases %>%
        fortify() %>%
        mutate(total = C1 + C2,
               total.inf = C1.q2.5 + C2.q2.5,
               total.sup = C1.q97.5 + C2.q97.5) %>%
        ggplot(aes(x = Index)) +
        geom_line(aes(y = C1), col = color1, lty =2) +
        geom_line(aes(y=C2), col = color2, lty =2) +
        geom_line(aes(y=total)) +
        geom_ribbon(aes(ymin = C1.q2.5, ymax = C1.q97.5), fill = "blue", alpha = 0.3) +
        geom_ribbon(aes(ymin = C2.q2.5, ymax = C2.q97.5), fill = "red", alpha = 0.3) +
        geom_ribbon(aes(ymin = total.inf, ymax = total.sup), fill = "black", alpha = 0.3) +
        geom_point(aes(y=covid.zoo)) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("New hospitalized cases") +
        xlab("")
    p2 <- sol$freq %>%
        fortify() %>%
        ggplot(aes(t_plot)) +
        geom_point(aes(y = freq)) +
        geom_errorbar(aes(ymin=obs_freq_lower, ymax=obs_freq_upper), width = 2) +
        geom_line(aes(y = pred.freq)) +
        geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.3) +
        scale_x_date(date_labels = "%b-%Y") +
        xlim(tmin, tmax) +
        theme_bw() +
        ylab("P.1 relative frequency") +
        xlab("date")
    grid.arrange(p1, p2, ncol = 1, nrow = 2)
}


## Ad hoc helper functions to convert back and forth parameters to logit scale ##
## Convert guesses that are probabilities to logit scale
prep.guess <- function(guess){
    if ("prevalence" %in% names(guess))
        guess["prevalence"] <- logit(guess["prevalence"]/100)
    guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- logit(guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(guess)
}
## convert guess back to prob. scale
prep.guess.inv <- function(guess) {
    if ("prevalence" %in% names(guess))
        guess["prevalence"] <- 100 * inv.logit(guess["prevalence"])
    guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- inv.logit(guess[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(guess)
}

## Convert coefficients that are at the logit scale back to
## proportions
prep.coef <- function(fit){
    fit.cf <- coef(fit)
    fit.parms <- fit.cf
    if ("prevalence" %in% names(fit.parms))
        fit.parms["prevalence"] <- inv.logit(fit.parms["prevalence"])*100
    if("IHR.V2.PROP" %in% names(fit.parms))
        fit.parms["IHR.V2.PROP"] <- exp(fit.parms["IHR.V2.PROP"])
    fit.parms[c("INIT.V2.FRAC", "PROB.REINFEC.V2")] <- inv.logit(fit.parms[c("INIT.V2.FRAC", "PROB.REINFEC.V2")])
    return(fit.parms)
}
## Convert confidence intervals that are at the logit scale
## back to proportions
prep.conf <- function(fit, ...){
    fit.conf <- confint(fit, ...)
    fit.conf[c("INIT.V2.FRAC", "PROB.REINFEC.V2"),] <- inv.logit(fit.conf[c("INIT.V2.FRAC", "PROB.REINFEC.V2"),])
    if ("prevalence" %in% names(fit.conf))
        fit.conf["prevalence",] <- inv.logit(fit.conf["prevalence",])*100
    return(fit.conf)
}
