### Overview

# This R-script provides code demonstrating how to fit the Cox proportional
# hazards model to a synthetic dataset simulated to resemble the Moderna
# mRNA-1273 data. This model includes several covariates - minority status, high
# risk status, and a risk score described in Gilbert, et al. (2021) - in
# addition to a time-varying predicted log10 antibody titre. The rate of decay
# for antibody titre is estimated from a separate seat of measurements from the
# Doria-Rose assay, as described in the supplementary material. The synthetic 
# datasets used for demonstration purposes in this script are simulated based on 
# summary statistics estimated from the mRNA-1273 and Doria-Rose datasets. 

### Setup
# load packages
library(brms)
library(dplyr)
library(mgcv)
library(survival)
library(tidyr)

## functions to simulate data 
# antibody decay data
simulate_abx_data = function() {
    dat = data.frame(ptid = rep(1:34, each = 3),
                     day = c(57, 119, 209),
                     sex = c(rep("M", 15 * 3), rep("F", 19 * 3)),
                     age = rep(rnorm(34, 52, 20), each = 3),
                     log_titre = c(rnorm(15, c(2.55, 2.32, 1.90), c(0.35, 0.42, 0.46)),
                                   rnorm(19, c(2.41, 2.08, 1.77), c(0.19, 0.26, 0.33))))
    return(dat)
}

# vaccine trial data 
simulate_cove_data = function() {
    dat = data.frame(ptid = 1:1615,
                     stratum = unlist(mapply(rep, times = c(22,20,24,23,18,26,441,171,143,188,168,146,189,36), x = 1:14)),
                     trt = unlist(mapply(rep, times = c(574, 1041), x = 0:1)),
                     time = pmax(round(unlist(c(rgamma(574 - 441, 8.3, 1/10), rgamma(441, 4.6, 1/10), rgamma(1041 - 36, 9.1, 1/10), rgamma(36, 5, 1/10)))), 1),
                     event = unlist(mapply(rep, times = c(574 - 441, 441, 1041 - 36, 36), x = rep(0:1, 2))),
                     d57titre = c(rep(0.08, 574), c(rnorm(1041 - 36, 2.4, 0.42), rnorm(36, 2.2, 0.42))), # LoD for placebo
                     minority = c(sample(c(rep(0,59),rep(1,74)), 574-441),
                                     sample(c(rep(0,319), rep(1,122)), 441),
                                     sample(c(rep(0,466), rep(1,539)), 1041 - 36),
                                     sample(c(rep(0,26), rep(1,10)), 36)),
                     highrisk = c(sample(c(rep(0,78),rep(1,55)), 574-441),
                                  sample(c(rep(0,340), rep(1,101)), 441),
                                  sample(c(rep(0,612), rep(1,393)), 1041 - 36),
                                  sample(c(rep(0,23), rep(1,13)), 36)),
                     riskscore = unlist(mapply(rnorm, 
                                               n = c(22,20,24,23,18,26,441,171,143,188,168,146,189,36),
                                               mean = c(-1.2,-0.2,0.1,-0.9,0.8,0.4,0.3,-1.4,0.0,0.0,-0.9,0.8,0.4,0.4),
                                               sd = c(1.1,0.9,1.1,0.7,0.9,0.6,0.9,1.0,1.4,1.1,0.5,0.8,0.5,0.8))),
                     weight = unlist(mapply(rep,times = c(22,20,24,23,18,26,441,171,143,188,168,146,189,36),
                                                x = c(22.4,36.3,107.3,128.1,65.4,199.3,3.1,5.4,14.9,18.2,9.0,29.4,1.5,1.3))))
    
    return(dat)
}

# conversion to start-stop format - use slope estimated from Doria-Rose data
# in example, everyone enrolled at day 0, not so in COVE. We work on calendar time
make_dat_long = function(dat, slope = -0.0043) {
    dat %>% 
        mutate(daystart = 0) %>% 
        group_by(ptid) %>% 
        reframe(tstart = seq(daystart, time - 1, by = 1),
                tstop = tstart + 1,
                event = case_when(tstop == time ~ event, TRUE ~ 0),
                trt = unique(trt),
                stratum = unique(stratum),
                minority = unique(minority),
                highrisk = unique(highrisk),
                riskscore = unique(riskscore),
                weight = unique(weight),
                day57titre = unique(d57titre),
                log_titre = unique(day57titre) + slope * (tstart - daystart)) %>% 
        ungroup()
}

### Antibody decay model

# simulate a dataset
abx_sim = simulate_abx_data()
cove_sim = make_dat_long(simulate_cove_data(), slope = -0.0043)

# fit antibody decay model
abx_fit = abx_sim %>% 
    mutate(l_cens = ifelse(log_titre == 1, -1, 0), # for measurements below the LoD
           day = day / 100) %>% # for numerical stability - remember to multiply x100 to get slope in days
    brm(formula = bf(log_titre | cens(l_cens) ~ day + sex + age + (1 | ptid)),
        family = gaussian(),
        data = .,
        prior = c(set_prior("normal(0, 2.5)", coef = "age", class = "b"),
                  set_prior("normal(0, 2.5)", coef = "sexM", class = "b"),
                  set_prior("normal(0, 2.5)", coef = "day", class = "b"),
                  set_prior("exponential(1)", class = "sd"),
                  set_prior("exponential(1)", class = "sigma")),
        chains = 4, 
        iter = 5e3)

# a non-Bayesian version of the Tobit model, fit via mgcv, for those interested
abx_fit2 = gam(log_titre ~ day + age + sex + s(ptid, bs = "re"), data = abx_sim, family = cnorm())


# fit Cox model to simulated COVE data
fit_cox = 
    coxph(Surv(tstart, tstop, event) ~ 
          minority + highrisk + riskscore + trt + trt:log_titre,
      data = cove_sim,
      id = ptid,
      robust = TRUE,
      weights = cove_sim$weight,
      model = TRUE)

# calculate VE at sequence of titres
VE_seq = 1 - exp(coef(fit_cox)["trt"] + coef(fit_cox)["trt:log_titre"] * seq(1.6,3.6,by=0.1))
plot(seq(1.6,3.6,by=0.1), VE_seq, "l")
