#!/usr/bin/env Rscript

### RI-v5.R
### last edited: 11 Jan 2021
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
ri.data <- read.csv("../../../../data/Rhode_Island/RI_formatted_20210606-TT.csv")
ri.data <- data.frame(ri.data)

## removing days before 61
idx.remove <- which(ri.data$daynum < 61)
if(length(idx.remove) > 0){
  ri.data <- ri.data[-idx.remove, ]
}

n.days <- nrow(ri.data)
days <- ri.data$daynum

### mcmc initialization (don't change)
end.day <- max(days, na.rm = TRUE) + 1
num.days <- end.day - 60 ## number of days with data

source("../../data.process.R")
dp <- data.process(ri.data)
# str(dp)

### Create spline expansions

# Cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis=round(num.days/7))

# create basis evaluation function matrices
Z <- eval.basis(bspl, 61:end.day)

# # cubic spline with multiple rates
# bspl.rr <- create.bspline.basis(c(61,end.day),
#                                 #nbasis=,
#                                 breaks=c(61, 84, 88, 92, 100, 130, 160, 190, 330, end.day),
#                                 norder=4)
# piecewise spline (starting at 61, ending at end.day) with split at day 100
# bspl.rr <- create.bspline.basis(rangeval = c(61, end.day),
#                                 breaks = c(61, 92, 122, 336, 398, end.day), norder = 1)
# Z.rr <- eval.basis(bspl.rr,61:end.day)

### create I-spline basis functions, evaluated every day from 61:end.day
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
# rr.start <- c(0.084, 0.310, 0.462, 0.545, 0.690, 0.7, 0.8, 0.85, 0.9, 0.9, 0.8, 0.8) 
# rr.start <- c(0.2, 0.5, 0.7, 0.8, 0.9)
# rr.start <- c(0.084, 0.310, 0.262, 0.145, 0.090) 
# rr.start <- c(0.30822, 0.2442, 0.09108, 0.07026, 0.04508, 0.0104, 0.1097)
rr.start <- c(0.30722, 0.22024, 0.177406, 0.10516, 0.052549, 0.023474, 0.035296)
rr.daily <- Z.rr %*% rr.start
# plot(rr.daily)

# starting values for beta (length should == ncol(Z))
# beta.strt <- c(1.07, 0.24, 0.10, 0.11, 0.07, 0.09, 0.03,
#                0.05, 0.03, 0.02, 0.03, 0.04, 0.01, 0.05,
#                0.06, 0.09, 0.10, 0.09, 0.09, 0.11, 0.10,
#                0.11, 0.06, 0.11, 0.06, 0.09, 0.07, 0.14,
#                0.07, 0.11, 0.08, 0.10, 0.10, 0.11, 0.13,
#                0.15, 0.15, 0.16, 0.14, 0.10, 0.10, 0.12,
#                0.16, 0.14, 0.15, 0.12, 0.14, 0.12, 0.16,
#                0.20, 0.30, 0.29, 0.29, 0.23, 0.19, 0.20,
#                0.25, 0.20, 0.25, 0.21, 0.21, 0.21, 0.26,
#                0.23, 0.19, 0.20)
beta.strt <- c(1.151, 0.451, 0.077, 0.088, 0.104, 
               0.058, 0.03, 0.054, 0.028, 0.019, 
               0.033, 0.027, 0.015, 0.03, 0.03, 
               0.062, 0.06, 0.094, 0.119, 0.086, 
               0.133, 0.073, 0.103, 0.057, 0.118,
               0.056, 0.073, 0.124, 0.118, 0.11,
               0.154, 0.072, 0.153, 0.121, 0.129, 
               0.115, 0.144, 0.109, 0.123, 0.084, 
               0.132, 0.067, 0.149, 0.089, 0.192, 
               0.031, 0.103, 0.147, 0.096, 0.1, 
               0.121, 0.196, 0.218, 0.19, 0.178, 
               0.222, 0.201, 0.181, 0.202, 0.136, 
               0.12, 0.118, 0.19, 0.096, 0.204,
               0.206)
# plot(Z %*% beta.strt)
# starting values for odesim parameters and other parameters
# prms.start <- c(10.20723, 0.01229, 0.04987, 0.20022, 
#                 0.64257, 
#                 0.67002, 0.81751, 
#                 0.39589, 0.28857, 0.36916, 
#                 155.29976, 328.03921, 
#                 0.01221, 0.02227, 0.03173, 0.04887, 0.07821, 0.14376, 0.2863, 0.30709, 
#                 2.34159, 2.9288, 2.93702, 2.43651, 1.83924, 1.10224, 1.07107, 2.17684, 
#                 2.90864, 2.44081, 1.593, 1.0191, 0.74885, 0.42019, 0.33322, 0.51494, 
#                 3.22507, 2.59056, 1.4468, 0.83452, 0.65603, 0.37747, 0.27688, 0.33386, 
#                 159.90631, 315.50196)
prms.start <- c(10.595277, 0.010079, 0.048401, 0.215648, 
                0.711116, 
                0.66643, 0.938992, 
                0.524219, 0.255902, 0.344169, 
                155.278723, 346.958336,

		0.010894, 0.018796, 0.034260, 0.053050, 0.089632, 0.165127, 0.311067, 0.334274, 
                
                2.466489, 3.098611, 2.957575, 2.452915, 1.854954, 1.12589, 1.092449, 2.177055, 
                2.861411, 2.337288, 1.650131, 1.082244, 0.864771, 0.493546, 0.369244, 0.538514, 
                3.708274, 3.057589, 1.318687, 0.660672, 0.502561, 0.28356, 0.199465, 0.151454, 
                159.144923, 379.495235)

names(prms.start) <- c("mean-time-vent", "death-prob-home-60", "death-prob-home-70","death-prob-home-80",
                       "tv-dev-len-hospstay",
                       "prob-icu-vent", "dev-ventdeath-mid",
                       "tv-dev-icu-frac_1","tv-dev-icu-frac_2", "tv-dev-icu-frac_3",
                       "tv-dev-icu-frac-endday_1", "tv-dev-icu-frac-endday_2",

			"tv-hosp-frac-10","tv-hosp-frac-20","tv-hosp-frac-30","tv-hosp-frac-40","tv-hosp-frac-50","tv-hosp-frac-60","tv-hosp-frac-70", "tv-hosp-frac-80",

                       "tv-contact-rate-10_1", "tv-contact-rate-20_1", "tv-contact-rate-30_1", "tv-contact-rate-40_1", "tv-contact-rate-50_1", "tv-contact-rate-60_1", "tv-contact-rate-70_1", "tv-contact-rate-80_1",
                       "tv-contact-rate-10_2", "tv-contact-rate-20_2", "tv-contact-rate-30_2", "tv-contact-rate-40_2", "tv-contact-rate-50_2", "tv-contact-rate-60_2", "tv-contact-rate-70_2", "tv-contact-rate-80_2",
                       "tv-contact-rate-10_3", "tv-contact-rate-20_3", "tv-contact-rate-30_3", "tv-contact-rate-40_3", "tv-contact-rate-50_3", "tv-contact-rate-60_3", "tv-contact-rate-70_3", "tv-contact-rate-80_3",
                       "tv-contact-rate-endday_1", "tv-contact-rate-endday_2")

## prior limits for other parameters (all uniform priors)
params.prior.min <- c(7, 0.001, 0.01, 0.1, 
                      0.1, 
                      0.4, 0.4, 
                      0.01, 0.01, 0.01,
                      130, 180,

			.001, .001, .005, .005, .01, .05, .1, .1,

                      .1, .1, .1, .1, .1, .1, .1, .1,
                      .1, .1, .1, .1, .1, .1, .1, .1,
                      .1, .1, .1, .1, .1, .1, .1, .1,
                      120, 250)

params.prior.max <- c(14, 0.2, 0.3, 0.4, 
                      2.0, 
                      1.0, 1.5,
                      1.5, 2.0, 2.0,
                      165, 360,

			.1, .1, .15, .15, .2, .2, .4, .4,

                      10, 10, 10, 10, 10, 10, 10, 10,
                      10, 10, 10, 10, 10, 10, 10, 10,
                      10, 10, 10, 10, 10, 10, 10, 10,
                      190, 410)


## names of parameters that are NOT odesim parameters, but are still to be estimated
non.odesim.params <- NULL
## inputs to odesim that are to be FIXED, and not estimated

# const.prms.start <- c("", 2)
# names(const.prms.start) <- c("symp-frac-davies", "steps-per-day")

## fixed hosp-frac: arithmetic mean of combined posterior from RI 0911, CT 0916, MA 0920 runs with different hosp-fracs
const.prms.start <- c("", 2,
                      4,
                      "347 354 361 368 375 382 389 396 403 410 417 424 431 438 445 452 459 466 473 480 487 494 501 508 515 522 529",
                      "0 0 0 0 1 3 10 44 50 55 55 68 53 73 92 69 177 421 369 699 1199 1864 5001 1949 2087 6411 6683 0",
                      "0 0 0 2 489 867 543 1477 1219 1456 1509 1252 1137 1604 1603 1078 2969 4500 2336 3573 5319 5466 13418 5175 4060 2380 2522 0",
                      "0 0 0 4 888 1279 513 1969 1497 1615 1402 1563 1291 2086 1895 1594 3998 4931 2871 3810 5149 5123 11424 4437 3332 1933 1959 0",
                      "0 0 0 3 700 1071 498 2072 1568 1680 1453 1751 1660 2312 1864 1928 4937 5732 3268 3891 5497 6851 10225 4565 2879 1799 1669 0",
                      "0 0 0 4 822 1368 717 2778 2060 2403 1882 2286 2236 3241 3445 3077 7381 10687 5349 6773 12405 9778 8805 3790 2825 1632 1440 0",
                      "0 0 0 3 565 987 646 2184 1610 2166 1916 2243 2505 4284 15636 13599 14099 11012 8089 8142 7015 3226 3280 1663 1474 866 776 0",
                      "0 0 0 1 71 142 294 873 630 1096 1385 2198 3135 5793 18837 14759 9358 3834 979 1013 1155 827 826 511 473 271 236 0",
                      "0 0 0 0 6 24 581 1140 817 1393 1378 2786 3332 3671 9535 6270 2473 1268 364 307 381 260 325 217 190 119 109 0")
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
prms.start >= params.prior.min & prms.start <= params.prior.max

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

## starting values for NB dispersion parameters
## first = new case counts
## second = new hospitalization counts
## third = total deaths (if lik.tot.deaths==TRUE)
## third = hosp deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## fourth = home deaths (if lik.hosp.deaths==TRUE and lik.home.deaths==TRUE)
## ...so nb.disp.start with have different lengths depending on what data are used
nb.disp.start <- c(1, 1, 0.8, 0.01, 1.1)



# create a "list" of starting values for all chains
betas.0 <- list()
ode.params.0 <- list()
const.params.0 <- const.prms.start
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
#                    const.params=const.params.0,
#                    non.odesim.params=non.odesim.params,
#                    s2.hosp.start=s2.hosp.0[[1]],
#                    s2.vent.start=s2.vent.0[[1]],
#                    s2.icu.start=s2.icu.0[[1]],
#                    lik.params.start=lik.params.0[[1]],
#                    Sigma.tune=Sigma.tune.0,
#                    var.tune=var.tune.0,
#                    ## mcmc.odesim options below
#                    df=ri.data, ## data
#                    loc="RI",
#                    lik.tot=FALSE, ## only use total new case data, not age-structured
#                    lik.age=TRUE,
#                    lik.hosp.new=TRUE, ## use daily new hospitalization data
#                    lik.curr = TRUE,
#                    lik.hosp.curr=TRUE,
#                    lik.icu.curr=TRUE,
#                    lik.vent.curr=TRUE,
#                    lik.tot.deaths=FALSE,
#                    lik.age.deaths=TRUE,
#                    lik.home.deaths = TRUE,
#                    lik.hosp.deaths = TRUE,
#                    lik.hosp.discharges=TRUE,
#                    total.size.constraint = TRUE, ## forces sympt to be close to data
#                    fixed.nb.disp = FALSE,
#                    active.surv = FALSE,
#                    p.asympt = 0.4, ## change if active.serv if TRUE...
#                    spline.beta=bspl, ## fda spline object for time-varying betas
#                    spline.rr=Z.rr, ## fda spline object for time-varying betas
#                    ode.params.prior.min=params.prior.min, ## prior minima  odesim inputs
#                    ode.params.prior.max=params.prior.max, ## prior maxima  odesim inputs
#                    start.day=61,
#                    end.day=max(ri.data$daynum,na.rm=TRUE)+1,
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
fit <- multichain.mcmc.odesim(parallel.type = "psock", ## "psock" for aci or "doMC"
                              n.chains = n.chs, ## number of independent MCMC chains
                              n.cores = n.chs, ## number of cores to use (keep equal to n.chains)
                              n.iter = 800, ## number of times to repeat n.mcmc iterations
                              n.mcmc.per.iter = 1000, ## number of mcmc iterations to run before saving
                              adapt.iter = 100, ## number of iterations between tuning (100 or 1000 is good)
                              thin.rt = 10, ## only save every "thin.rt"-th iteration
                              resample = FALSE, ## ignore for now
                              df = ri.data, ## data
                              loc = "RI",
                              save.file.name = paste(output_dir,"out",sep="/"), #"~/scratch/covid19/20210111/RI-jan6/out",  ## output file directory
                              inf.file.dir = "../../",
                              odesim.ver = "v5", ## either "v4" or "v5"
                              introday = 55, ## day of start of epidemic
                              betas.lst = betas.0, ## starting values
                              ode.params.lst = ode.params.0, ## starting values
                              lik.params.lst = lik.params.0, ## starting values
                              rr.params.lst = rr.params.0, ## starting values
                              s2.hosp.lst = s2.hosp.0, ## starting values
                              s2.vent.lst = s2.vent.0, ## starting values
                              s2.icu.lst = s2.icu.0, ## starting values
                              Sigma.tune = Sigma.tune.0,
                              var.tune = var.tune.0, ## starting value for proposal variance.  Can adjust
                              lik.tot = FALSE, ## only use total new case data, not age-structured
                              lik.age = TRUE, ## use age-structured new case data
                              lik.hosp.new = TRUE, ## use daily new hospitalization data
                              lik.curr = TRUE,
                              lik.hosp.curr = TRUE, ## use current hospitalization data
                              lik.icu.curr = TRUE, ## no current ICU data
                              lik.vent.curr = TRUE, ## No current on-vent data
                              lik.tot.deaths = FALSE, ## total deaths
                              lik.age.deaths = TRUE, ## age-structured deaths
                              lik.home.deaths = TRUE,
                              lik.hosp.deaths = TRUE,
                              lik.hosp.discharges = TRUE, ## use new hosp discharges
                              total.size.constraint = TRUE, ## forces sympt to be close to data
                              fixed.nb.disp = FALSE,
                              active.surv = FALSE,
                              p.asympt = 0.4, ## change if active.serv if TRUE...
                              spline.beta = bspl, ## fda spline object for time-varying betas
                              spline.rr = Z.rr, ## fda spline object for time-varying reporting rate
                              ode.params.prior.min = params.prior.min, ## prior minima  odesim inputs
                              ode.params.prior.max = params.prior.max, ## prior maxima  odesim inputs
                              const.params = const.params.0,
                              non.odesim.params = non.odesim.params,
                              start.day = 61, ## don't change this!
                              end.day = max(ri.data$daynum,na.rm=TRUE) + 1, ## timing info
                              odepath = odepath, ## path to odesim
                              p.vecs = delay.probs, ## testing delay probabilities
                              adapt = "ShabyWells", ## adaptive tuning type.  "None" means no tuning.
                              plot.save = FALSE, ## write out plots occasionally
                              sf.choice = FALSE,
                              mult = 1
)


