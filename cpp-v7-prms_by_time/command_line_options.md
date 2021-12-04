### File describing the command-line usage of the executable 'odesim'

##### First argument is always the output filename

You can call `./odesim myoutput.tdl`

or `./odesim none`

and the 'none' keyword will skip outputting to any file, and everything gets sent to stdout instead.

##### -binary-output

With the `-binary-output` flag `odesim` will out output ODE results as
raw doubles, rather than tab-delimited text.  This is significantly faster but
requires the user reshape the data based on the number of state variables.
e.g.
```python
import subprocess
import numpy as np

NCOLS = 307  # 34*9 compartments + time

proc = subprocess.run(['./odesim', 'none', '-binary-output'],
                      stdout=subprocess.PIPE)

arr = np.frombuffer(proc.stdout).reshape((-1,307))
``` 

##### -beta

You can call the function with an arbitrary number of beta parameters (transmission parameters) that allow for time-varying transmission in the model.  For example, you can call

   `./odesim none -beta 1.0 0.9 0.8 0.7`
   
and this will run the model with four different transmission rates, evenly distributed, in the order above, for the period between March 1 and the end of the simulation which is denoted by `tf = time final`.  


> WARNING
>
> For now, you have to set time final on the command line before you assign the betas
>
> e.g. `./odesim none -tf 120 -beta 3 4 5 6`



##### -contact-rate-10, -contact-rate-20, -contact-rate-30, -contact-rate-40, -contact-rate-50, -contact-rate-60, -contact-rate-70, -contact-rate-80 
These command-line option were **deprecated** on Jan 22 2021. See **-tv-contact-rate-x0** options
    
This is the relative contact rate (or population mixing level) by age group. `-contact-rate-10` is the rate at which individuals in the 10-19 age group are contacting or mixing with other individuals, irrespective of the age of the contacts.  All mixing or contact rates are relative to the rate for 0-9 individuals.

```diff
! Recommended prior distributions 

!  -contact-rate-10 = [0.3 - 3.0]
!  -contact-rate-20 = [0.3 - 3.0]
!  -contact-rate-30 = [0.3 - 3.0]
!  -contact-rate-40 = [0.3 - 3.0]
!  -contact-rate-50 = [0.3 - 3.0]
!  -contact-rate-60 = [0.3 - 3.0]
!  -contact-rate-70 = [0.3 - 3.0]
!  -contact-rate-80 = [0.3 - 3.0]

```

##### -contact-rate-postld-10, -contact-rate-postld-20, -contact-rate-postld-30, -contact-rate-postld-40, -contact-rate-postld-50, -contact-rate-postld-60, -contact-rate-postld-70, -contact-rate-postld-80
These command-line option were **deprecated** on Jan 22 2021. See **-tv-contact-rate-x0** options

As above, but these are the population's "new" contact rates after the lockdown

```diff
! Recommended prior distributions 

!  -contact-rate-postld-10 = [0.3 - 3.0]
!  -contact-rate-postld-20 = [0.3 - 3.0]
!  -contact-rate-postld-30 = [0.3 - 3.0]
!  -contact-rate-postld-40 = [0.3 - 3.0]
!  -contact-rate-postld-50 = [0.3 - 3.0]
!  -contact-rate-postld-60 = [0.3 - 3.0]
!  -contact-rate-postld-70 = [0.3 - 3.0]
!  -contact-rate-postld-80 = [0.3 - 3.0]

```




##### -death-prob-home-60, -death-prob-home-70, -death-prob-home-80

This is the probability of dying at home for individuals in their 60s, 70s, and 80s.

```diff
! Recommended prior distributions 

!  -death-prob-home-60 = [0.0 - 0.10]
!  -death-prob-home-70 = [0.0 - 0.20]
!  -death-prob-home-80 = [0.0 - 0.30]

```

##### -death-prob-nonicu-80

This is the probaility that a hospitalized (non-ICU) individual in the 80+ age group dies as a result of infection.  Default is 0.05.  If you change this, then then half of this probability is used for the 70-79 age group.

```diff
! Recommended prior distributions 

!  -death-prob-nonicu-80 = [0.02 - 0.30]

```


##### -death-prob-postvent

This is the probaility that an extubated individual dies while still in the ICU.  Default is 0.189 which is set from a single study.  This probability is distributed among the 60+ age groups.

```diff
! Recommended prior distributions 

!  -death-prob-postvent = [0.0 - 0.30]

```



##### -death-rate-output

This will skip all ODE integration, and using the currently set parameters will ouput the following 30 numbers in this order: the infection fatality rate (IFR) for the nine age groups, the symptomatic case fatality rate (sCFR) for the nine age groups, the hospitalization fatality rate (HFR) for the nine age groups, the age-adjusted IFR for the entire population, the age-adjusted sCFR for the entire population, and the age-adjusted HFR for the entire population.


##### -dev-icu-frac, -dev-icu-frac-phase2
These command-line option were **deprecated** on Feb 17 2021. See **-tv-dev-icu-frac** options

These are the deviations in the parameters that determines the fraction of hospitalized individuals that are admitted to the ICU. This deviation will apply across all age groups.  Defaults are 1.0.  If you set this higher than 2.0, it will be reset back down to 2.0.  The phase2 deviation is a new deviation that applies to a new epidemic phase with better clinical management (fewer patients going to ICU).  The start date of phase2 is set below.  If these two parameters are the same value, this means that the probability of progression to the ICU does not change from phase 1 to phase 2.

```diff
! Recommended prior distribution [0.5 - 1.5]
```

##### -dev-icu-frac-phase2beginday

Begin day of improvement clinical management during the COVID-19 epidemic phase, in principle, the period during summer 2020 when dexamethasone and prone positioning began to be used to reduce the likelihood of ICU admission.  Default is set to 153 which is June 1 2020.



##### -dev-len-hospstay
This command-line option was **deprecated** on Sep 06 2021. See **-tv-dev-len-hospstay** options

This is the deviation in the parameter that determines the length of the typical hospital stay on the medical floor.  Default is 1.0.  If you set this higher than 2.5, it will be reset back down to 2.5 

```diff
! Recommended prior distribution [0.5 - 2.0]
```


##### -dev-ventdeath-mid

This is the deviation in the parameter that determines the fraction of 40-70 individuals that die after being on a ventilator. Default is 0.7.  Do not set this higher than 1.7. 

```diff
! Recommended prior distribution [0.5 - 1.3]
```


##### -dev-ventdeath-70, -dev-ventdeath-80 

This is the deviation in the parameter that determines the fraction of 70-79 and 80+ individuals that die after being on a ventilator. Defaults are 1.0.  

```diff
! Strict bounds for prior distributions (do not make wider)

!  -dev-ventdeath-70 = [0.8 - 1.4]
!  -dev-ventdeath-80 = [0.7 - 1.1]

```



##### -diag

This is a flag that prints out some basic diagnostic info as to who was infected, who was hospitalized, and who died.  To be used like this:

    `./odesim none -diag`
    
    
    
    
    
##### -firstlockdown-endday
These command-line option were **deprecated** on Jan 22 2021. See **-tv-contact-rate-beginday** and **-tv-contact-rate-endday** options

This is the last day of the "full part" of the first lockdown, which is set to 121 as a default which is April 30.  The -contact-rate-postld-nn parameters change on this day.



##### -hosp-frac-10, -hosp-frac-20, -hosp-frac-30, -hosp-frac-40, -hosp-frac-50, -hosp-frac-60, -hosp-frac-70, -hosp-frac-80
These command-line option were **deprecated** on Mar 08 2021. See **-tv-hosp-frac-x0** options
    
This is the fraction of symptomatic individuals that are eventually hospitalized, by age group. `-hosp-frac-10` is the hospitalization probability for individuals in the 10-19 age group. Once this probability is set, it is also applied to the 0-9 age group. `-hosp-frac-20` is for the 20-29 age group, `-hosp-frac-30` is for the 30-39 age group, etc.

```diff
! Recommended prior distributions 

!  -hosp-frac-10 = [0.0 - 0.2]
!  -hosp-frac-20 = [0.0 - 0.2]
!  -hosp-frac-30 = [0.0 - 0.2]
!  -hosp-frac-40 = [0.0 - 0.5]
!  -hosp-frac-50 = [0.0 - 0.5]
!  -hosp-frac-60 = [0.0 - 1.0]
!  -hosp-frac-70 = [0.0 - 1.0]
!  -hosp-frac-80 = [0.0 - 1.0]

```


    
##### -introday

The introduction day of the first case. Feb 1 2020 is day 32, March 1 2020 is day 61, etc.


##### -loc

The location currently being modeled. E.g.

    `./odesim none -loc PA`
    
for Pennsylvania. Acceptable arguments currently are PA, MA, and RI.  Default is RI.


##### -mean-time-vent

This is the mean time that a surviving patient spends on a ventilator.  Default is set to 10.8 days.


##### -min-mixinglevel-00, -min-mixinglevel-10, -min-mixinglevel-20, -min-mixinglevel-30, -min-mixinglevel-40, -min-mixinglevel-50, -min-mixinglevel-60, -min-mixinglevel-70, -min-mixinglevel-80

This is the minimum allowable mixing level per age group. Default is for these to be set to zero. Argument must be between 0.0 and 1.0. For example, if you set

   `./odesim none -min-mixinglevel-80 0.30`

this means that individuals in the 80+ age group will only be allowed to have their mixing level reduced by 70% during a lockdown, even if other age groups have their mixing reduced by 80% or 90% or more.


##### -prob-icu-vent

This is the probability that a patient in the ICU will progress to mechanical ventilation.  Reasonable values are 0.75 to 0.99.


##### -rel-beta-hosp

This is the relative beta, i.e. the relative population mixing parameter, for a hospitalized individual (including ICU and ventilated). Default is 0.2. Must be between 0.0 and 1.0. Note that this beta does not change when social distancing is enforced.

```diff
! Recommended prior distribution [0.0 - 0.3]
```


##### -remove99

This option (or flag) removes the first 99 columns of the output.  These are the columns that keep track of the susceptible individuals, exposed individuals, and asymptomatic individuals.  These individuals are never enterred into the likelihood function, so it is not necessary to output them.



##### -self-isolation-factor

A factor describing the extent to which symptomatic individuals self-isolate by not leaving home and contacting others.  1.0 means that they behave identically to non-symptomatic individuals.  0.0 means that they self-isolate immediately after the first symptom appears and do not contact anyone until either (1) they need to be hospitalized or (2) all symptoms have disappeared.


##### -steps-per-day

This is the number of steps the Runge-Kutta integrator takes per day.  Default is set to 100.  Can set to 5 or 10 to improve speed in trial inference runs.



##### -susc-0-20, -susc-60-100

This is the relative susceptibility of certain age groups. Normally, the 30-39 age group is our reference group.  Here, the entire 20-59 age group is the reference group, since there is unlikely to be enough information to differentiate susceptibility among 10-year age bands.  The susceptibility of the 20-59 age group is 1.0, and the two other age-group susceptibilities may be able to be inferred.  Default settings are that the susceptibility of the 0-20 age group is 0.6 (slight adjustment from an odds-value of 0.51) and the susceptibility of the 60+ age group is 1.0.

```diff
! Recommended prior distributions 

!  -susc-0-20   = [0.2 - 1.5]
!  -susc-60-100 = [0.5 - 1.5]

```



##### -symp-frac

This command-line option was deprecated on July 8 2020.



##### -symp-frac-equal

This is the probability that infections progress to symptoms. With this option it is set to be equal for all ages, e.g.

    `./odesim none -symp-frac-equal 0.601`

is set as the probability for all age classes, where 0.601 is the current (2020/07/08) simple weighted average across all studies we have reviewed that distinguished between asymptomatic and pre-symptomatic patients.


##### -symp-frac-davies

This command-line option is a flag and takes no other arguments, e.g.

    `./odesim none -symp-frac-davies`

and this flag sets the probabilities (by age) of progressing to symptoms as the posterior means inferred in Davies et al (Nature Medicine, June 16 2020).




##### -tf

Time at which the ODEs are stopped. Day 1 is Jan 1 2020.  So, if you want to run the simulation through to April 30 2020, you would call

   `./odesim none -tf 121`
   
The 'time initial' right now is set to zero, which means that it is set to Jan 1 2020.  Cases are introduced via the `-introday` command-line argument.  


##### -time-symp-to-hosp

This is the average time (in days) between symptoms appearance and admission to hospital. It is the same for all ages right now.

```diff
! Recommended prior distribution [2.0 - 10.0]
```



##### -tv-contact-rate-10, -tv-contact-rate-20, -tv-contact-rate-30, -tv-contact-rate-40, -tv-contact-rate-50, -tv-contact-rate-60, -tv-contact-rate-70, -tv-contact-rate-80

Each of these options take a space-separated list of non-negative floating-point numbers. All lists must have identical length. This list implementation of contact rates is meant to replace the old `-contact-rate-x0` and `-contact-rate-postld-x0` options in version 5. Each set of 8 contact rates (1 from each age group) indicates a phase in contact patterns.

This is the relative contact rate (or population mixing level) by age group. `-tv-contact-rate-10` is the rate at which individuals in the 10-19 age group are contacting or mixing with other individuals, irrespective of the age of the contacts.  All mixing or contact rates are relative to the rate for 0-9 individuals (fixed at to 1.0).

##### -tv-contact-rate-beginday, -tv-contact-rate-endday

The begining and ending day for each phase in contact patterns (as defined by `-tv-contact-rate-x0` options). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.




##### -tv-contact_coeff-00, -tv-contact_coeff-10, -tv-contact_coeff-20, -tv-contact_coeff-30, -tv-contact_coeff-40, -tv-contact_coeff-50, -tv-contact_coeff-60, -tv-contact_coeff-70, -tv-contact_coeff-80

Each of these options take a space-separated list of non-negative floating-point numbers. All lists must have identical length. This list implementation of contact coefficients is meant to be used with [Belgium CoMix](https://doi.org/10.1038/s41598-020-78540-7) contact matrices. The first set of 9 contact coefficients (1 from each age group) will be used with matrix from CoMix wave 1 survey. The other subsequent sets will be used with matrix from CoMix wave 5 survey.

`-tv-contact-rate-*` and `-tv-contact_coeff-*` options are mutually exclusive.

##### -tv-contact_coeff-beginday, -tv-contact_coeff-endday

The begining and ending day for each phase in contact coefficients (as defined by `-tv-contact_coeff-x0` options). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.


##### -tv-dev-icu-frac
This option takes a space-separated list of non-negative floating-point numbers. This implementation is meant to replace the old `-dev-icu-frac` and `-dev-icu-frac-phase2` options in version 5 and to allow more than two values for `dev-icu-frac`. 

##### -tv-dev-icu-frac-beginday, -tv-dev-icu-frac-endday

The begining and ending day for each phase in `dev-icu-frac` patterns (as defined by `-tv-dev-icu-frac` option). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.



##### -tv-dev-len-hospstay
This option takes a space-separated list of non-negative floating-point numbers. This implementation is meant to replace the old `-dev-len-hospstay` options in version 5 and to allow more than one value for `dev-len-hospstay`. 


##### -tv-dev-len-hospstay-beginday, -tv-dev-len-hospstay-endday

The begining and ending day for each phase in `dev-len-hospstay` patterns (as defined by `-tv-dev-len-hospstay` option). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.


##### -tv-hosp-frac-10, -tv-hosp-frac-20, -tv-hosp-frac-30, -tv-hosp-frac-40, -tv-hosp-frac-50, -tv-hosp-frac-60, -tv-hosp-frac-70, -tv-hosp-frac-80

Each of these options take a space-separated list of non-negative floating-point numbers. All lists must have identical length. This list implementation of hospitalization fractions is meant to replace the old `-hosp-frac-x0` options in version 5. Each set of 8 hosp-frac (1 from each age group) indicates a phase in hospitalization patterns.

##### -tv-hosp-frac-beginday, -tv-hosp-frac-endday

The begining and ending day for each phase in hospitalization patterns (as defined by `-tv-hosp-frac-x0` options). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.



##### -tv-vaccinees-10, -tv-vaccinees-20, -tv-vaccinees-30, -tv-vaccinees-40, -tv-vaccinees-50, -tv-vaccinees-60, -tv-vaccinees-70, -tv-vaccinees-80

Each of these options take a space-separated list of non-negative floating-point numbers. All lists must have identical length. Each set of 8 numbers of fully vaccinated individuals (1 from each age group) indicates a phase in the vaccination campaign. We assume that days within the same phase of the campaign have identical pace of vaccination (i.e. number of vaccinees per day). The number of vaccinees aged 0-9 is 0 throughout the simulation.

##### -tv-vaccinees-beginday, -tv-vaccinees-endday

The begining and ending day for each phase in the vaccination campaign (as defined by `-tv-vaccinees-x0` options). Each of these options take a space-separated list of non-negative increasing integer numbers. All lists should have identical length. If not, new elements will be appended to the shorter list according to the corresponding elements in the longer list. 

Each begining-ending day pair must satisfy the relationship: begining day < ending day.

#### -no-vaccine-wastage

Flag to turn off vaccine wastage, i.e. allowing vaccinating seropositive individuals who might not known they already contracted the virus.

#### -rr 

Average reporting rate (to be used in determining which fraction of the R-class might be getting the vaccine).
