#!/usr/bin/env Rscript

### MA-v5.R
### last edited: 09 Dec 2020
### (edited from Dr. Hank's example script)

args = commandArgs(trailingOnly=TRUE) 
if (length(args)==0) {
  output_dir = "output"
} else {
  output_dir = args[1]
}

# Number of mcmc chains to run in parallel (Note: MUST BE LESS THAN 20!)
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
## for iSpline
library(splines2)

### compile the odesim model
# odepath <- "../../../../cpp-v5-discharges-nonhospdeaths/"
odepath="../../../../cpp-v7-prms_by_time/"

## uncomment to recompile odesim
# system(paste("rm ",odepath,"odesim",sep=""))
# system(paste("make -C", odepath))

### read in data
# ma.data <- read.csv("../../../../data/Massachusetts/Past_data/MA_20200906_data_7daysmoothed.csv")
ma.data <- read.csv("../../../../data/Massachusetts/MA_20210606_data_7daysmoothed-HI_TT_full.csv")
# ma.data <- read.csv("../../../../data/Massachusetts/MA_20210606_data-TT.csv")

source("../../../code/data.process.R")

## removing days before 61
idx.remove <- which(ma.data$daynum < 61)
if(length(idx.remove) > 0){
  ma.data <- ma.data[-idx.remove,]
}


n.days <- nrow(ma.data)
days <- ma.data$daynum


dp <- data.process(ma.data, loc="MA")


### mcmc initialization (don't change)
end.day <- max(days, na.rm = TRUE) + 1
num.days <- end.day - 60 ## number of days with data

### Create spline expansions

# Cubic spline with one basis function every 7 days
bspl <- create.bspline.basis(c(61, end.day), nbasis = round(num.days/7))
## zero-order spline with one rate (constant reporting rate)
# bspl.rr <- create.bspline.basis(c(61,end.day), nbasis = 1, norder = 1)
## cubic spline with multiple rates
# bspl.rr <- create.bspline.basis(c(61,end.day),
#                                 #nbasis=,
#                                 breaks=c(61, 84, 92, 122, 155, 190, 280, 336, 398, end.day),
#                                 norder=4)
## piecewise spline (starting at 61, ending at end.day) with split at day 100
# bspl.rr <- create.bspline.basis(rangeval = c(61, end.day),
#                                 breaks = c(61, 92, 122, 336, 398, end.day), norder = 1)
# create basis evaluation function matrices
Z <- eval.basis(bspl, 61:end.day)

# Z.rr <- eval.basis(bspl.rr, 61:end.day)

## i-spline
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
delay.probs <- list(p1 = c(1,rep(0,7)),
                    p2 = c(1,rep(0,7)),
                    p3 = c(1,rep(0,7)),
                    p4 = c(1,rep(0,7)))


###############################################################################
### 2. Starting values for parameters
###############################################################################

## starting values for reporting rates
# rr.start <- c(0.184, 0.310, 0.462, 0.545, 0.590, 0.6, 0.7, 0.7, 0.8, 0.85, 0.85, 0.9)
# rr.start <- c(0.2, 0.5, 0.7, 0.8, 0.9)
# rr.start <- c(0.30822, 0.1442, 0.07108, 0.01026, 0.00508, 0.0064, 0.21997) # c(0.3, 0.05, 0.05, 0.15, 0.09, 0.1, 0.14)
rr.start <- c(0.169356, 0.319222, 0.165428, 0.012537, 0.016388, 0.009203, 0.197095)

# rr.start <- 0.75 
rr.daily <- Z.rr %*% rr.start
# plot(rr.daily)

# starting values for beta (length should == ncol(Z))
# beta.strt <- c(1.85743, 0.57785, 0.26205, 0.08928, 0.0365, 
#                0.08816, 0.03625, 0.03177, 0.01695, 0.03711, 
#                0.01852, 0.0312, 0.02767, 0.01857, 0.21653, 
#                0.13382, 0.55594, 0.241, 0.29191, 0.30139, 
#                0.35535, 0.25238, 0.24789, 0.39559, 0.27738, 
#                0.31338, 0.34441, 0.35395, 0.41954, 0.44227, 
#                0.47667, 0.31325, 0.30038, 0.43502, 0.55595, 
#                0.63475, 0.62583, 0.58574, 0.60448, 0.66829, 
#                0.41069, 0.55621, 0.75292, 0.91601, 0.50287, 
#                0.48744, 0.61693, 0.3133, 0.15832, 0.55478,
#                0.56205, 0.62839, 0.75493, 0.75204, 0.60609, 
#                0.73222, 0.63526, 0.6982, 0.51794, 0.56078, 
#                0.30154, 0.19266, 0.10432, 0.08174, 0.12143, 
#                0.13938)
beta.strt <- c(1.8088, 0.5822, 0.1935, 0.0855, 0.0288, 
               0.1016, 0.0206, 0.0372, 0.012, 0.0374, 
               0.0132, 0.031, 0.0253, 0.0313, 0.2316, 
               0.1797, 0.5793, 0.2552, 0.2851, 0.3896,
               0.3894, 0.2691, 0.2976, 0.4232, 0.3435,
               0.2773, 0.4552, 0.3599, 0.4154, 0.5309, 
               0.4697, 0.3842, 0.3483, 0.513, 0.5548, 
               0.6185, 0.6858, 0.5913, 0.5926, 0.6316, 
               0.5608, 0.2879, 0.4225, 0.9265, 0.895, 
               0.8979, 0.6279, 0.3537, 0.3056, 0.6301,
               0.6579, 0.7102, 0.7899, 0.5989, 0.856, 
               0.7144, 0.3503, 0.8227, 0.7258, 0.815, 
               0.5256, 0.2651, 0.0696, 0.0738, 0.206, 
               0.2718)
# plot(61:end.day, Z %*% beta.strt)


## starting values for odesim parameters and other parameters
# prms <- c(10.99667, 0.01086, 0.0503, 0.28931, 
#           0.86869, 0.84738, 
#           262.00422, 
#           1.15648, 0.35954, 0.28552, 
#           147.49649, 254.28894, 
#           0.79474, 0.77854, 
#           0.01184, 0.0192, 0.03882, 0.05856, 0.10146, 0.18602, 0.34059, 0.36572, 
#           1.37801, 2.71566, 2.68002, 2.29713, 1.80621, 1.33558, 1.50108, 3.78051, 
#           1.16911, 1.20597, 0.83848, 0.59371, 0.38292, 0.22142, 0.18327, 0.27677, 
#           0.95866, 1.02364, 0.70951, 0.52308, 0.37553, 0.2191, 0.17816, 0.2341, 
#           1.11015, 1.02076, 0.62826, 0.40568, 0.28377, 0.1549, 0.11837, 0.1721, 
#           149.48321, 299.65458, 393.36795)
prms <- c(11.017294, 0.010379, 0.047463, 0.292739, 
          0.80305,
          1.244906, 0.390731, 0.305416, 
          147.83073, 256.7339, 
          0.790264, 0.781007,

	0.010894, 0.018796, 0.034260, 0.053050, 0.089632, 0.165127, 0.311067, 0.334274, 
          
          1.493209, 2.821223, 2.782355, 2.392119, 1.875654, 1.394608, 1.535098, 3.858751, 
          1.063458, 1.149609, 0.800348, 0.568948, 0.367774, 0.213771, 0.174833, 0.261816, 
          0.888936, 0.985566, 0.687154, 0.50911, 0.366675, 0.21564, 0.173765, 0.227295, 
          1.037401, 0.973343, 0.605049, 0.394203, 0.276852, 0.152704, 0.115944, 0.166566, 
          149.399613, 299.60515, 391.920998)

names(prms) <- c("mean-time-vent", "death-prob-home-60", "death-prob-home-70", "death-prob-home-80", 
                 "tv-dev-len-hospstay",
                 "tv-dev-icu-frac_1","tv-dev-icu-frac_2","tv-dev-icu-frac_3",
                 "tv-dev-icu-frac-endday_1","tv-dev-icu-frac-endday_2",
                 "prob-icu-vent", "dev-ventdeath-mid", 
                 
		"tv-hosp-frac-10","tv-hosp-frac-20","tv-hosp-frac-30","tv-hosp-frac-40","tv-hosp-frac-50","tv-hosp-frac-60","tv-hosp-frac-70", "tv-hosp-frac-80",

                 "tv-contact-rate-10_1","tv-contact-rate-20_1","tv-contact-rate-30_1","tv-contact-rate-40_1","tv-contact-rate-50_1","tv-contact-rate-60_1","tv-contact-rate-70_1","tv-contact-rate-80_1",
                 "tv-contact-rate-10_2","tv-contact-rate-20_2","tv-contact-rate-30_2","tv-contact-rate-40_2","tv-contact-rate-50_2","tv-contact-rate-60_2","tv-contact-rate-70_2","tv-contact-rate-80_2",
                 "tv-contact-rate-10_3","tv-contact-rate-20_3","tv-contact-rate-30_3","tv-contact-rate-40_3","tv-contact-rate-50_3","tv-contact-rate-60_3","tv-contact-rate-70_3","tv-contact-rate-80_3",
                 "tv-contact-rate-10_4","tv-contact-rate-20_4","tv-contact-rate-30_4","tv-contact-rate-40_4","tv-contact-rate-50_4","tv-contact-rate-60_4","tv-contact-rate-70_4","tv-contact-rate-80_4",
                 "tv-contact-rate-endday_1", "tv-contact-rate-endday_2", "tv-contact-rate-endday_3")


## prior limits for other parameters (all uniform priors)
params.prior.min=c(7.0, 0.001, 0.01, 0.1, 
                   0.1,
                   0.01, 0.01, 0.01,
                   110, 250,
                   0.4, 0.4,
		
		.001, .001, .005, .005, .01, .05, .1, .1,
                   
                   .1, .1, .1, .1, .1, .1, .1, .1,
                   .1, .1, .1, .1, .1, .1, .1, .1,
                   .1, .1, .1, .1, .1, .1, .1, .1,
                   .1, .1, .1, .1, .1, .1, .1, .1,
                   100, 230, 330)

params.prior.max=c(14.0, 0.2, 0.3, 0.4, 
                   2.0,
                   2.0, 2.0, 2.0,
                   200, 380,
                   1.0, 1.5, 

		.1, .1, .15, .15, .2, .2, .4, .4,
                   
                   10, 10, 10, 10, 10, 10, 10, 10,
                   10, 10, 10, 10, 10, 10, 10, 10,
                   10, 10, 10, 10, 10, 10, 10, 10,
                   10, 10, 10, 10, 10, 10, 10, 10,
                   190, 305, 420)

## rough calculation for variance
range=params.prior.max-params.prior.min
ode.params.var.start=range^2/(mean(range^2))
names(ode.params.var.start)=names(prms)

## starting values for non ode and constant parameters
non.odesim.params = NULL

## fixed hosp-frac: arithmetic mean of combined posterior from RI 0911, CT 0916, MA 0920 runs with different hosp-fracs
const.prms.start=c("", 2,
                   4,
                   "364 371 378 385 392 399 406 413 420 427 434 441 448 455 462 469 476 483 490 497 504 511 518 525 532",
                   "0 1 73 50 260 316 594 761 1245 933 1067 1000 1400 1300 2200 5500 7000 12900 15400 23000 25200 25000 19000 65400 56800 0",
                   "0 114 5089 3414 5890 8104 11260 13313 14618 12837 11361 9700 7100 8800 17600 34900 32800 45300 47600 64200 51800 51900 38000 34200 20000 0",
                   "0 262.30 7295.87 4262.15 7001.69 10440.36 12856.47 16457.12 17931.13 14963.25 14600.10 13794.38 11861.08 16354.71 29626.57 49377.62 39136.34 49273.12 49847.89 65732.33 47966.83 47548.82 33963.45 30096.84 16615.96 0",
                   "0 239.70 6667.13 3894.85 6398.31 9540.64 11748.53 15038.88 16385.87 13673.75 13341.90 12605.62 10838.92 14945.29 27073.43 45122.38 35763.66 45026.88 45552.11 60067.67 43833.17 43451.18 31036.55 27503.16 15184.04 0",
                   "0 302.36 6421.90 3455.52 7674.71 11654.15 16508.68 17385.39 20844.83 16785.03 20586.38 22143.13 22950.78 31633.04 45026.61 67573.56 65285.21 74640.52 65419.82 69727.30 39171.13 33517.56 25710.26 21537.39 11980.17 0",
                   "0 262.64 5578.10 3001.48 6666.29 10122.85 7400.79 13304.32 17623.59 15180.12 27705.30 36443.52 48511.61 66421.82 64171.44 61272.81 52890.35 51597.17 40961.56 41090.43 23805.31 20793.52 18115.33 13965.33 7942.69 0",
                   "0 18.13 377.76 283.61 3830.21 4627.25 7712.89 7061.97 13743.49 29800.61 45812.12 45011.99 42377.18 49766.45 41376.92 30098.11 21665.27 16083.03 9896.04 8095.06 5516.23 5216.96 6073.19 4099.10 2491.16 0",
                   "0 12.87 268.24 201.39 2719.79 3285.75 9027.64 5792.32 14899.09 49401.24 63159.20 46701.36 19560.43 9678.69 8325.03 10355.52 6159.17 5279.28 3722.57 3587.21 2707.33 2571.96 2301.23 2098.18 1285.98 0")
names(const.prms.start) = c("symp-frac-davies", "steps-per-day", 
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


############################################################
##
## Empirical start values for NB dispersion parameters
## (added 8/21/20)
############################################################


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
hosp.nb.size.hat=1/median(abs((resids.hosp^2-hosp.hat))/hosp.hat^2)
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
nb.disp.start=c(sympt.nb.size.hat,hosp.nb.size.hat,deaths.nb.size.hat)
nb.disp.start



# ## checking visual fit
# source("../../traj.process.R")
# source("../../traj.from.params.R")
# source("../../plots.odesim.R")
# tp=traj.process(traj.from.params(beta=beta.strt,params=prms,const.params=const.prms.start,non.odesim.params=non.odesim.params,
#                                  loc="MA",odepath=odepath,tf=end.day),loc="MA",odesim.version="v5")
# par(mfrow=c(2,3))
# plots.odesim(dp,tp,rr=rr.daily[-1])
# plot(tp$tot.hosp.new.odesim)

# create a "list" of starting values for all chains
betas.0=list()
ode.params.0=list()
rr.params.0=list()
lik.params.0=list()
s2.hosp.0=list()
s2.vent.0=list()
s2.icu.0=list()
for(k in 1:n.chs){
  betas.0[[k]]=beta.strt
  ode.params.0[[k]]=prms
  rr.params.0[[k]]= rr.start # rep(.7,ncol(Z.rr))
  lik.params.0[[k]]=nb.disp.start
  s2.hosp.0[[k]]=1
  s2.vent.0[[k]]=1
  s2.icu.0[[k]]=1
}


##
## proposal variance object
##
## order= 
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
#                    s2.hosp.start=s2.hosp.0[[1]],
#                    s2.vent.start=s2.vent.0[[1]],
#                    s2.icu.start=s2.icu.0[[1]],
#                    lik.params.start=lik.params.0[[1]],
#                    Sigma.tune=Sigma.tune.0,
#                    var.tune=var.tune.0,
#                    ## mcmc.odesim options below
#                    df=ma.data, ## data
#                    loc="MA",
#                    lik.tot = FALSE, ## only use total new case data, not age-structured
#                    lik.age = TRUE, ## use age-structured new case data
#                    lik.curr = TRUE,
#                    lik.hosp.new = TRUE, ## use daily new hospitalization data
#                    lik.hosp.curr = TRUE, ## use current hospitalization data
#                    lik.icu.curr = TRUE, ## use current ICU data
#                    lik.vent.curr = TRUE, ## use current on-vent data
#                    lik.tot.deaths = FALSE, ## use total deaths
#                    lik.home.deaths = FALSE,
#                    lik.hosp.deaths = FALSE,
#                    lik.age.deaths = TRUE,
#                    total.size.constraint = TRUE,
#                    fixed.nb.disp = FALSE,
#                    spline.beta=bspl, ## fda spline object for time-varying betas
#                    spline.rr=Z.rr, #bspl.rr, ## fda spline object for time-varying betas
#                    ode.params.prior.min=params.prior.min, ## prior minima  odesim inputs
#                    ode.params.prior.max=params.prior.max, ## prior maxima  odesim inputs
#                    start.day=61,
#                    end.day=end.day,
#                    introday=55, ## timing info
#                    odepath=odepath, ## path to odesim
#                    odesim.ver="v5",
#                    p.vecs=delay.probs, ## testing delay probabilities
#                    thin=10, ## thinning rate for saving posterior samples
#                    plot.rate=100,
#                    print.iter=FALSE,
#                    c0=1,
#                    c1=.8,
#                    adapt.type="ShabyWells",
#                    adapt.iter=100,
#                    # hosp.report.rate=0.4,
#                    const.params=const.prms.start,
#                    non.odesim.params=non.odesim.params,
#                    sf.choice=FALSE
# )


###############################################################################
### 3. MCMC chains in parallel
###############################################################################

n.iter = 800
n.mcmc.per.iter = 1000

# library(doMC)
# registerDoMC(cores=n.chs)

# source("../../multichain.mcmc.odesim.2.0.R")
source("../../multichain.mcmc.odesim.2.0-v7-TT_CT.R")

### Run MCMC
fit <- multichain.mcmc.odesim(parallel.type = "psock", ## "psock" for aci or "doMC"
                              n.chains = n.chs, ## number of independent MCMC chains
                              n.cores = n.chs, ## number of cores to use (keep equal to n.chains)
                              n.iter = n.iter, ## number of times to repeat n.mcmc iterations
                              n.mcmc.per.iter = n.mcmc.per.iter, # number of mcmc iterations to run before saving
                              resample = FALSE, ## ignore for now
                              df = ma.data, ## data
                              loc = "MA",
                              save.file.name = paste(output_dir,"out",sep="/"), #"~/scratch/covid19/20210111/MA-jan6/out",  ## output file directory
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
                              lik.curr = TRUE,
                              lik.hosp.new = TRUE, ## use daily new hospitalization data
                              lik.hosp.curr = TRUE, ## use current hospitalization data
                              lik.icu.curr = TRUE, ## use current ICU data
                              lik.vent.curr = TRUE, ## use current on-vent data
                              lik.tot.deaths = FALSE, ## use total deaths
                              lik.home.deaths = FALSE,
                              lik.hosp.deaths = FALSE,
                              lik.age.deaths = TRUE,
                              lik.hosp.discharges=FALSE,
                              total.size.constraint = TRUE,
                              fixed.nb.disp = FALSE,
                              spline.beta = bspl, ## fda spline object for time-varying betas
                              spline.rr = Z.rr, #bspl.rr, ## fda spline object for time-varying reporting rate
                              ode.params.prior.min = params.prior.min, ## prior minima  odesim inputs
                              ode.params.prior.max = params.prior.max, ## prior maxima  odesim inputs
                              start.day = 61, ## don't change this!
                              end.day = max(ma.data$daynum,na.rm=TRUE)+1, ## timing info
                              odepath = odepath, ## path to odesim
                              p.vecs = delay.probs, ## testing delay probabilities
                              adapt = "ShabyWells", ## adaptive tuning type.  "None" means no tuning.
                              adapt.iter = 100, ## number of iterations between tuning (100 or 1000 is good)
                              thin.rt = 10, ## only save every "thin.rt"-th iteration
                              plot.save = FALSE, ## write out plots occasionally,
                              const.params = const.prms.start,
                              non.odesim.params = non.odesim.params,
                              sf.choice = FALSE
)
