
#!/usr/bin/env Rscript

### CT-v5.R
### last edited: 23 Jan 2021
### Ephraim Hanks, Nathan Wikle, Emily Strong

args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0) {
  output_dir = "output"
} else {
  output_dir = args[1]
}

### Number of mcmc chains to run in parallel (Note: MUST BE LESS THAN 20!)
n.chs = 5

###############################################################################
### 1. Preliminaries
###############################################################################

### necessary packages

## library for vectorized multinomial
library(mc2d)
## library for multivariate normal distribution
library(mvtnorm)
## library to read in .xlsx files (not needed to read in .csv files)
library(readxl)
## library for spline-based expansions
library(fda)
## packages for parallelization:
library(snow); library(doParallel); library(foreach)
## package for multivariate normal
library(msm)
## package for iSplines:
library(splines2)

### compile the odesim model
# odepath="../../../../cpp-v5-discharges-nonhospdeaths/"
odepath="../../../../cpp-v7-prms_by_time/"

## uncomment to recompile odesim
# system(paste("rm ",odepath,"odesim",sep=""))
# system(paste("make -C", odepath))

### read in data
ct.data <- read.csv("../../../../data/Connecticut/Connecticut_covid_data_with_tests_CDC_20210606-TT.csv")

ct.data <- data.frame(ct.data)

n.days <- nrow(ct.data)
days <- ct.data$daynum

### mcmc initialization (don't change)
end.day <- max(days, na.rm = TRUE) + 1
num.days <- end.day - 60 ## number of days with data

# fix ct mixup (169:171 -> 49251, 49236, 49289)
# ct.data[169:171, 5]
# ct.data[169:171, 5] <- c(49236, 49251, 49289)
# 
# # some cumulative age sympt data are decreasing...
# ct.data[c(101:106, 108:109, 115:116, 134:135, 155:158,
#                 141:144, 169:170, 218:219), 16:24] <- NA
# 
# # some cumulative age death data are decreasing...
# ct.data[c(92:113, 134:137, 142, 148:149, 158, 166:169, 177:178,
#   191:201, 213:215, 219, 233, 239, 255:256, 277:278, 303), 34:42] <- NA
# 
# # remove last four days of death data
# ct.data[309:312, 7:8] <- NA
# ct.data[309:312, 34:42] <- NA

# change column to match other data
colnames(ct.data)[43] <- "Cumulative_hospital_deaths"
colnames(ct.data)[44] <- "Cumulative_hospital_discharges"

source("../../data.process.R")
dp <- data.process(ct.data)
str(dp)
  
# No vent data!
# No new hosp data for 0-10 and 10-20
# Removed last four days of death data (incomplete)
# No home or hosp death data
# No (new) hospital discharge data 
#   (we do have cumulative data, although it looks incomplete...)


### Create spline expansions

# Cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days / 7))
# create basis evaluation function matrices
Z <- eval.basis(bspl, 61:end.day)


# # cubic spline with multiple rates
# bspl.rr <- create.bspline.basis(c(61,end.day),
#                                 #nbasis=,
#                                 breaks=c(61, 84, 92, 100, 130, 160, 190, 330, 339, end.day),
#                                 norder=4)
# piecewise spline (starting at 61, ending at end.day) with split at day 100
# bspl.rr <- create.bspline.basis(rangeval = c(61, end.day),
#                                 breaks = c(61, 92, 122, 336, 398, end.day), norder = 1)

# Z.rr <- eval.basis(bspl.rr,61:end.day)

### create I-spline basis functions, evaluated every day from 61:end.day
### only i-spline through August, then piecewise-constant...
# knots <- c(84, 92, 100, 108, 130, 160, 190)
# Z.rr <- iSpline(x = 61:end.day, knots = knots, degree = 2)
# Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
# Z.rr <- Z.rr[,c(1,4:7)]
# Z.rr <- iSpline(x = 61:end.day, knots = c(84, 92, 122, 155, 190, 280, 336), degree = 2)
# Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
# Z.rr <- Z.rr[,c(1,3:8)]

Z.rr <- iSpline(x = 61:end.day, knots = c(92, 122, 183, 288, 336, 381, 410), degree = 2)
Z.rr <- cbind(rep(1, nrow(Z.rr)), Z.rr)
Z.rr <- Z.rr[,c(1,3:8)]


### Creating testing delay probabilities

# list of 4 vectors, each corresponding to a reporting period:
#   p1: Beginning - March 13,
#   p2: March 14 - March 23,
#   p3: March 24 - March 27,
#   p4: March 28 - present
# the 1st term = prob of 0 day delay, 2nd = prob. of 1 day delay, etc.
delay.probs <- list(p1 = c(1, rep(0,7)),
                    p2 = c(1, rep(0,7)),
                    p3 = c(1, rep(0,7)),
                    p4 = c(1, rep(0,7)))


###############################################################################
### 2. Starting values for parameters
###############################################################################

# starting values for reporting rate
### b-spline
# rr.start <- c(0.0402744, 0.066826, 0.7472345, 0.9714792, 0.9732023, 0.960018, 0.9742731, 0.9778607, 0.9948033, 0.995954, 0.9902845, 0.9802008)
# rr.start <- c(0.03889, 0.08, 0.4503, 0.6932, 0.72142, 0.74067, 0.7624, 0.79465, 0.859, 0.8632, 0.84, 0.84)
### piecewise
# rr.start <- c(0.2, 0.5, 0.77, 0.80, 0.85)
### i-spline
# rr.start <- c(0.084, 0.310, 0.262, 0.105, 0.089)
# rr.start <- c(0.30822, 0.1442, 0.07108, 0.01026, 0.00508, 0.0064, 0.21997)
# 0.004471	0.377583	0.548901	0.021194	0.011802	0.0093	0.007496

rr.start <- c(0.40822, 0.13042, 0.20108, 0.10026, 0.05508, 0.0264, 0.01297)
rr.daily <- Z.rr %*% rr.start


# starting values for beta (length should == ncol(Z))
# beta.strt <- c(0.55902, 0.11276, 0.02383, 0.02924, 0.05633, 
#                0.01592, 0.01394, 0.01048, 0.00658, 0.00878, 
#                0.00738, 0.00284, 0.01199, 0.01965, 0.07478, 
#                0.10171, 0.08407, 0.13069, 0.05384, 0.10349, 
#                0.09313, 0.07055, 0.07613, 0.12199, 0.12128, 
#                0.09983, 0.12811, 0.08222, 0.10812, 0.15562, 
#                0.12634, 0.13997, 0.15284, 0.12513, 0.17308, 
#                0.21021, 0.19301, 0.15981, 0.23151, 0.16721, 
#                0.14003, 0.16690, 0.16001, 0.2202, 0.24724, 
#                0.25104, 0.18053, 0.17556, 0.19748, 0.216, 
#                0.25364, 0.20343, 0.26277, 0.26797, 0.3066, 
#                0.22926, 0.26638, 0.25545, 0.18024, 0.17569, 
#                0.17039, 0.14213, 0.10904, 0.10783, 0.17593, 
#                0.3054)
beta.strt <- c(1.139, 0.02, 0.009, 0.019, 0.025, 
               0.019, 0.012, 0.008, 0.008, 0.006, 
               0.007, 0.002, 0.006, 0.032, 0.027, 
               0.181, 0.111, 0.073, 0.069, 0.103, 
               0.075, 0.103, 0.102, 0.083, 0.106, 
               0.138, 0.112, 0.1, 0.123, 0.185, 
               0.133, 0.155, 0.184, 0.106, 0.135, 
               0.215, 0.139, 0.172, 0.238, 0.163, 
               0.123, 0.162, 0.198, 0.186, 0.154, 
               0.154, 0.177, 0.134, 0.122, 0.221, 
               0.271, 0.235, 0.252, 0.316, 0.265, 
               0.316, 0.259, 0.276, 0.25, 0.219, 
               0.206, 0.234, 0.208, 0.253, 0.295, 
               0.335)
# plot(Z %*% beta.strt)

# starting values for odesim parameters and other parameters
# prms.start <- c(9.1015365, 0.0115819, 0.0505501, 0.2616012, 
#                 0.7270915, 
#                 0.5938406, 0.7392946, 
#                 0.8859089, 0.4648103, 0.5437361, 
#                 155.9033635, 309.0638943,
#                 2.3788176, 4.6975804, 5.1024163, 4.2984845, 3.7027998, 2.6327943, 2.5264236, 5.3710544, 
#                 2.7491518, 2.3457372, 1.5200176, 0.9835879, 0.6843779, 0.3605143, 0.3037911, 0.4525456, 
#                 2.3736681, 1.6066941, 1.1601296, 0.8193094, 0.5992292, 0.3293416, 0.2680385, 0.3726569, 
#                 143.728027, 331.5725625)

prms.start <- c(9.133498, 0.011578, 0.046256, 0.243863,
                0.803489, 
                0.623433, 0.744137, 
                0.524219, 0.395308, 0.566127, 
                155.738997, 311.856189,

		0.010894, 0.018796, 0.034260, 0.053050, 0.089632, 0.165127, 0.311067, 0.334274, 

                2.487262, 4.93022, 5.312914, 4.523388, 3.824248, 2.730508, 2.625011, 5.66409, 
                2.682156, 2.30017, 1.485833, 0.959471, 0.667819, 0.351964, 0.298454, 0.440301, 
                2.390019, 1.62525, 1.167718, 0.820438, 0.598332, 0.326729, 0.267877, 0.373643, 
                145.28292, 311.598776)

names(prms.start) <- c("mean-time-vent", "death-prob-home-60", "death-prob-home-70","death-prob-home-80",
                       "tv-dev-len-hospstay",
                       "prob-icu-vent", "dev-ventdeath-mid", 
                       
                       "tv-dev-icu-frac_1","tv-dev-icu-frac_2","tv-dev-icu-frac_3",
                       "tv-dev-icu-frac-endday_1","tv-dev-icu-frac-endday_2",
			
			"tv-hosp-frac-10","tv-hosp-frac-20","tv-hosp-frac-30","tv-hosp-frac-40","tv-hosp-frac-50","tv-hosp-frac-60","tv-hosp-frac-70", "tv-hosp-frac-80",

                       "tv-contact-rate-10_1","tv-contact-rate-20_1","tv-contact-rate-30_1","tv-contact-rate-40_1","tv-contact-rate-50_1","tv-contact-rate-60_1","tv-contact-rate-70_1","tv-contact-rate-80_1",
                       "tv-contact-rate-10_2","tv-contact-rate-20_2","tv-contact-rate-30_2","tv-contact-rate-40_2","tv-contact-rate-50_2","tv-contact-rate-60_2","tv-contact-rate-70_2","tv-contact-rate-80_2",
                       "tv-contact-rate-10_3","tv-contact-rate-20_3","tv-contact-rate-30_3","tv-contact-rate-40_3","tv-contact-rate-50_3","tv-contact-rate-60_3","tv-contact-rate-70_3","tv-contact-rate-80_3",
                       "tv-contact-rate-endday_1", "tv-contact-rate-endday_2")

## prior limits for other parameters (all uniform priors)
## restrict tv-dev-icu-frac_1 to the range of RI posterior
params.prior.min <- c(7, 0.001, 0.01, 0.1,  
                      0.1, 
                      0.4, 0.4,
                      0.4, 0.01, 0.01,
                      140, 280,

			.001, .001, .005, .005, .01, .05, .1, .1,

                      .1, .1, .1, .1, .1, .1, .1, .1,
                      .1, .1, .1, .1, .1, .1, .1, .1,
                      .1, .1, .1, .1, .1, .1, .1, .1,
                      100, 250)

params.prior.max <- c(14, 0.2, 0.3, 0.4, 
                      2.0, 
                      1.0, 1.5,
                      0.6, 2.0, 2.0,
                      210, 360,

			.1, .1, .15, .15, .2, .2, .4, .4, 
 
                      10, 10, 10, 10, 10, 10, 10, 10,
                      10, 10, 10, 10, 10, 10, 10, 10,
                      10, 10, 10, 10, 10, 10, 10, 10,
                      190, 390)


## names of parameters that are NOT odesim parameters, but are still to be estimated
non.odesim.params <- NULL
## inputs to odesim that are to be FIXED, and not estimated

## fixed hosp-frac: arithmetic mean of combined posterior from RI 0911, CT 0916, MA 0920 runs with different hosp-fracs
const.prms.start <- c("", 2,
                      4,
                      "406 413 420 426 433 440 447 456 463 470 477 484 491 498 505 512 519 526",
                      "0 4304 926 669 686 651 684 792 2635 2632 2632 5644 12167 17837 13677 9458 5024 33082 27449 11104 10538 5253 4702 5808 0",
                      "0 19755 3512 2231 2687 2372 2428 3219 10272 7980 7980 14357 28196 40311 29063 21232 11545 10886 11283 6910 9814 3846 3516 4199 0",
                      "0 28797 4905 2966 4215 3349 3292 5289 16082 10139 10139 16916 28356 38148 25144 19311 10715 8973 9638 6095 9011 3190 3016 3507 0",
                      "0 29747 5407 3239 4827 3812 3968 7530 18372 30897 9903 35480 27371 29590 20270 15769 9050 7539 7905 5027 8457 2637 2508 2938 0",
                      "0 37405 6970 4326 8996 10342 6897 21072 56865 77782 178 46206 24567 19382 14537 12486 7599 6633 7122 4379 7156 2365 2307 2586 0",
                      "0 29651 6889 4710 15597 34897 39027 39746 63357 25208.5 25208.5 20437 13053 10132 7950 7414 4476 4079 5083 2598 3510 1426 1483 1608 0",
                      "0 31655 19860 17097 17391 28728 34146 24130 20559 5228.5 5228.5 5348 3961 3340 2877 2917 1700 1707 2354 1001 1354 557 671 708 0",
                      "0 32665 25007 22147 13729 9467 8500 5575 5076 1724.5 1724.5 1956 1689 1469 1304 1371 762 895 1166 485 626 283 364 338 0")

names(const.prms.start) <- c("symp-frac-davies", "steps-per-day",
                             "time-symp-to-hosp",
                             "tv-vaccinees-endday",
                             "tv-vaccinees-10",
                             "tv-vaccinees-20",
                             "tv-vaccinees-30",
                             "tv-vaccinees-40",
                             "tv-vaccinees-50",
                             "tv-vaccinees-60",
                             "tv-vaccinees-70",
                             "tv-vaccinees-80")

## check to make sure priors and starting values line up
rbind(params.prior.min, prms.start, params.prior.max)

##
##
## checking visual fit (good for choosing prms and beta.strt above)
##
##
source("../../traj.process.R")
source("../../data.process.R")
source("../../traj.from.params-v7-TT.R")
# source("../../traj.from.params.R")
source("../../plots.odesim.R")

## rough calculation for tuning variance - make proposal variance smaller for parameters with smaller prior range
range <- params.prior.max-params.prior.min
ode.params.var.start <- range^2/(mean(range^2))
names(ode.params.var.start) <- names(prms.start)


# par(mfrow=c(2,2))
## sympt
sympt.loess = loess(dp$tot.sympt.new~dp$days)
sympt.hat=predict(sympt.loess)
resids.sympt=dp$tot.sympt.new[!is.na(dp$tot.sympt.new)]-sympt.hat
days.sympt=dp$days[!is.na(dp$tot.sympt.new)]
## estimating size parameter
sympt.nb.size.hat=1/median((resids.sympt^2-sympt.hat)/sympt.hat^2)
## check with plots and simulations
mean.sympt=sympt.hat
mean.sympt[mean.sympt<=0]=.0001
# plot(days.sympt,rnbinom(length(days.sympt),mu=mean.sympt,size=sympt.nb.size.hat),type="p",main="sympt")
# points(days.sympt,sympt.hat,type="l")
# points(dp$days,dp$tot.sympt.new,type="b",col="red")

## hosp
hosp.loess = loess(dp$tot.hosp.new~dp$days)
hosp.hat=predict(hosp.loess)
resids.hosp=dp$tot.hosp.new[!is.na(dp$tot.hosp.new)]-hosp.hat
days.hosp=dp$days[!is.na(dp$tot.hosp.new)]
## estimating size parameter
hosp.nb.size.hat=1/median((resids.hosp^2-hosp.hat)/hosp.hat^2)
## check with plots and simulations
mean.hosp=hosp.hat
mean.hosp[mean.hosp<=0]=.0001
# plot(days.hosp,rnbinom(length(days.hosp),mu=mean.hosp,size=hosp.nb.size.hat),type="p",main="hosp")
# points(days.hosp,hosp.hat,type="l")
# points(dp$days,dp$tot.hosp.new,type="b",col="red")

## deaths
deaths.loess = loess(dp$tot.deaths.new~dp$days)
deaths.hat=predict(deaths.loess)
resids.deaths=dp$tot.deaths.new[!is.na(dp$tot.deaths.new)]-deaths.hat
days.deaths=dp$days[!is.na(dp$tot.deaths.new)]
## estimating size parameter
deaths.nb.size.hat=1/median((resids.deaths^2-deaths.hat)/deaths.hat^2)
## check with plots and simulations
mean.deaths=deaths.hat
mean.deaths[mean.deaths<=0]=.0001
# plot(days.deaths,rnbinom(length(days.deaths),mu=mean.deaths,size=deaths.nb.size.hat),type="p",main="deaths")
# points(days.deaths,deaths.hat,type="l")
# points(dp$days,dp$tot.deaths.new,type="b",col="red")



## starting values for NB dispersion parameters
## first = new case counts
## second = new hospitalization counts
## third = total deaths (if lik.tot.deaths==TRUE)
## third = hosp deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fourth = home deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fifth = discharges
## ...so nb.disp.start with have different lengths depending on what data are used
nb.disp.start=c(sympt.nb.size.hat, hosp.nb.size.hat, deaths.nb.size.hat)
nb.disp.start

# nb.disp.start <- c(1, 1, 0.8) 

# create a "list" of starting values for all chains
betas.0 <- list()
ode.params.0 <- list()
rr.params.0 <- list()
lik.params.0 <- list()
s2.hosp.0 <- list()
s2.vent.0 <- list()
s2.icu.0 <- list()
for(k in 1:n.chs){
  betas.0[[k]] <- beta.strt
  ode.params.0[[k]] <- prms.start
  rr.params.0[[k]] <- rr.start
  lik.params.0[[k]] <- nb.disp.start
  s2.hosp.0[[k]] <- 1
  s2.vent.0[[k]] <- 1
  s2.icu.0[[k]] <- 1
}


# initialize proposal densities

### initialize proposal variance
var.tune.0=c(.0001,.001,.007,25,.01)
Sigma.tune.0=list(Sigma.tune.beta=diag(length(betas.0[[1]])),
                  Sigma.tune.ode.params=diag(ode.params.var.start),
                  Sigma.tune.rr.params=diag(length(rr.params.0[[1]])),
                  Sigma.tune.s2=diag(3),
                  Sigma.tune.ll.params=diag(length(lik.params.0[[1]]))
)

# OPT_USE_PY_LL=TRUE
# source("../../traj.from.params.R")
# source("../../loglik.odesim.4.0.R", chdir=TRUE)
# source("../../mcmc.odesim.2.0.R")
source("../../traj.from.params-v7-TT.R")
source("../../loglik.odesim.new.curr_CT.R", chdir=TRUE)
source("../../mcmc.odesim.2.0-v7-TT_CT.R")
source("../../traj.process.R")
source("../../data.process.R")
source("../../plots.odesim.R")
# fit <- mcmc.odesim(n.mcmc=601,
#                    beta.start=betas.0[[1]],
#                    report.rate.params.start=rr.params.0[[1]],
#                    ode.params.start=ode.params.0[[1]],
#                    const.params=const.prms.start,
#                    non.odesim.params=non.odesim.params,
#                    s2.hosp.start=s2.hosp.0[[1]],
#                    s2.vent.start=s2.vent.0[[1]],
#                    s2.icu.start=s2.icu.0[[1]],
#                    lik.params.start=lik.params.0[[1]],
#                    Sigma.tune=Sigma.tune.0,
#                    var.tune=var.tune.0,
#                    ## mcmc.odesim options below
#                    df=ct.data, ## data
#                    loc="CT",
#                    lik.old = FALSE,
#                    lik.tot=FALSE, ## only use total new case data, not age-structured
#                    lik.age=TRUE,
#                    lik.hosp.new=TRUE, ## use daily new hospitalization data
#                    lik.hosp.curr=TRUE,
#                    lik.icu.curr=TRUE,
#                    lik.vent.curr=FALSE,
#                    lik.curr=TRUE,
#                    lik.tot.deaths=FALSE,
#                    lik.age.deaths=TRUE,
#                    lik.home.deaths = FALSE,
#                    lik.hosp.deaths = FALSE,
#                    lik.hosp.discharges=FALSE,
#                    total.size.constraint = TRUE, ## forces sympt to be close to data
#                    fixed.nb.disp = FALSE,
#                    active.surv = FALSE,
#                    p.asympt = 0.4, ## change if active.serv if TRUE...
#                    spline.beta=bspl, ## fda spline object for time-varying betas
#                    spline.rr=Z.rr, ## fda spline object for time-varying betas
#                    ode.params.prior.min=params.prior.min, ## prior minima  odesim inputs
#                    ode.params.prior.max=params.prior.max, ## prior maxima  odesim inputs
#                    start.day=61,
#                    end.day=max(ct.data$daynum,na.rm=TRUE)+1,
#                    introday=55, ## timing info
#                    odepath=odepath, ## path to odesim
#                    odesim.ver="v5",
#                    p.vecs=delay.probs, ## testing delay probabilities
#                    thin=1, ## thining rate for saving posterior samples
#                    c0=1,
#                    c1=.8,
#                    adapt.type="ShabyWells",
#                    adapt.iter=100,
#                    s2.beta.start=.01,
#                    s2.rr.start=.01,
#                    t.adapt.start=0,
#                    prop.type="tnorm",
#                    plot.save=TRUE,
#                    plot.rate=200,
#                    print.iter=FALSE,
#                    plot.name="trace.plots.pdf",
#                    sf.choice = FALSE
# )
# str(fit)
## save(fit,file="fit.Rdata")



###############################################################################
### 3. MCMC chains in parallel
###############################################################################

# source("../../multichain.mcmc.odesim.2.0.R")
source("../../multichain.mcmc.odesim.2.0-v7-TT_CT.R")

# library(doMC)
# registerDoMC(cores=n.chs)

### Run MCMC
fit <- multichain.mcmc.odesim(n.chains = n.chs,
                              n.cores = n.chs,
                              n.iter = 800,
                              resample = FALSE,
                              parallel.type = "psock", ## "psock" for aci or "doMC"
                              inf.file.dir = "../../", ## location of mcmc.odesim.R / loglik.odesim.R / traj.from.params.R### remaining arguments go into mcmc.odesim
                              n.mcmc.per.iter = 1000, ## number of mcmc samples per iteration
                              save.file.name = paste(output_dir,"out",sep="/"), #"~/scratch/covid19/20210123/CT-jan6-n41/out",  ## output file directory,
                              betas.lst = betas.0, ## starting values
                              ode.params.lst = ode.params.0, ## starting values
                              lik.params.lst = lik.params.0, ## starting values
                              rr.params.lst = rr.params.0, ## starting values
                              s2.hosp.lst = s2.hosp.0, #var# starting values
                              s2.vent.lst = s2.vent.0, ## starting values
                              s2.icu.lst = s2.icu.0, ## starting values
                              Sigma.tune.start = Sigma.tune.0,
                              var.tune.start = var.tune.0,
                              adapt = "ShabyWells",
                              thin.rt = 10,
                              df = ct.data, ## data
                              loc = "CT",
                              odesim.ver = "v5", ## either "v4" or "v5"
                              introday = 55, ## day of start of epidemic
                              lik.old = FALSE,
                              lik.tot = FALSE, ## only use total new case data, not age-structured
                              lik.age = TRUE, ## use age-structured new case data
                              lik.hosp.new = TRUE, ## use daily new hospitalization data
                              lik.curr = TRUE,
                              lik.hosp.curr = TRUE, ## use current hospitalization data
                              lik.icu.curr = TRUE, ## no current ICU data
                              lik.vent.curr = FALSE, ## No current on-vent data
                              lik.tot.deaths = FALSE, ## total deaths
                              lik.age.deaths = TRUE, ## age-structured deaths
                              lik.home.deaths = FALSE,
                              lik.hosp.deaths = FALSE,
                              lik.hosp.discharges = FALSE, ## use new hosp discharges
                              total.size.constraint = TRUE, ## forces sympt to be close to data
                              fixed.nb.disp = FALSE,
                              active.surv = FALSE,
                              p.asympt = 0.4, ## change if active.serv if TRUE...
                              spline.beta = bspl, ## fda spline object for time-varying betas
                              spline.rr = Z.rr, ## fda spline object for time-varying reporting rate
                              ode.params.prior.min = params.prior.min, ## prior minima  odesim inputs
                              ode.params.prior.max = params.prior.max, ## prior maxima  odesim inputs
                              const.params = const.prms.start,
                              non.odesim.params = non.odesim.params,
                              start.day = 61, ## don't change this!
                              end.day = max(ct.data$daynum,na.rm=TRUE) + 1, ## timing info
                              odepath = odepath, ## path to odesim
                              p.vecs = delay.probs, ## testing delay probabilities
                              plot.save = FALSE, ## write out plots occasionally
                              sf.choice = FALSE,
                              mult = 1
)
