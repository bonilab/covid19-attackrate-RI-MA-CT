Location	MA	 
Last Day of Data Used	523	 
Day of MCMC Run	2021-11-27	 
 	 	 
 	 	 
Number of MCMC Chains	5	 
Number of Iterations Per Chain	592000	 
Burnin (values before this discarded)	20000	 
Adaptive Tuning Type	ShabyWells	 
lik.tot	FALSE	 
lik.age	TRUE	 
lik.hosp.new	TRUE	 
lik.hosp.curr	TRUE	 
lik.icu.curr	TRUE	 
lik.vent.curr	TRUE	 
lik.tot.deaths	FALSE	 
lik.home.deaths	FALSE	 
lik.age.deaths	TRUE	 
lik.hosp.deaths	FALSE	 
lik.hosp.discharges	FALSE	 
lik.curr	TRUE	 
lik.old	FALSE	 
 	 	 
odesim version	v5	 
 	 	 
Mean Log-Likelihood	-30905.0641483939	 
SD of Log-Likelihood	26.6174675698021	 
 	 	 
Dbar	61811.329762582	 
Dhat	62016.1667916765	 
DIC	61606.4927334874	 
DIC (Gelman)	63103.7362348951	 
 	 	 
 	 	 
Contact Rate Betas (spline length)	66	 
Contact Rate Betas (min day of spline)	61	 
Contact Rate Betas (max day of spline)	524	 
 	 	 
Reporting Rate (spline length)	7	 
Reporting Rate (min day of spline)	61	 
Reporting Rate (max day of spline)	524	 
 	 	 
 	 	 
Included ODESIM Params	Prior Min	Prior Max
mean-time-vent	7	14
death-prob-home-60	0.001	0.2
death-prob-home-70	0.01	0.3
death-prob-home-80	0.1	0.4
tv-dev-len-hospstay	0.1	2
tv-dev-icu-frac_1	0.01	2
tv-dev-icu-frac_2	0.01	2
tv-dev-icu-frac_3	0.01	2
tv-dev-icu-frac-endday_1	110	200
tv-dev-icu-frac-endday_2	250	380
prob-icu-vent	0.4	1
dev-ventdeath-mid	0.4	1.5
tv-hosp-frac-10	0.001	0.1
tv-hosp-frac-20	0.001	0.1
tv-hosp-frac-30	0.005	0.15
tv-hosp-frac-40	0.005	0.15
tv-hosp-frac-50	0.01	0.2
tv-hosp-frac-60	0.05	0.2
tv-hosp-frac-70	0.1	0.4
tv-hosp-frac-80	0.1	0.4
tv-contact-rate-10_1	0.1	10
tv-contact-rate-20_1	0.1	10
tv-contact-rate-30_1	0.1	10
tv-contact-rate-40_1	0.1	10
tv-contact-rate-50_1	0.1	10
tv-contact-rate-60_1	0.1	10
tv-contact-rate-70_1	0.1	10
tv-contact-rate-80_1	0.1	10
tv-contact-rate-10_2	0.1	10
tv-contact-rate-20_2	0.1	10
tv-contact-rate-30_2	0.1	10
tv-contact-rate-40_2	0.1	10
tv-contact-rate-50_2	0.1	10
tv-contact-rate-60_2	0.1	10
tv-contact-rate-70_2	0.1	10
tv-contact-rate-80_2	0.1	10
tv-contact-rate-10_3	0.1	10
tv-contact-rate-20_3	0.1	10
tv-contact-rate-30_3	0.1	10
tv-contact-rate-40_3	0.1	10
tv-contact-rate-50_3	0.1	10
tv-contact-rate-60_3	0.1	10
tv-contact-rate-70_3	0.1	10
tv-contact-rate-80_3	0.1	10
tv-contact-rate-10_4	0.1	10
tv-contact-rate-20_4	0.1	10
tv-contact-rate-30_4	0.1	10
tv-contact-rate-40_4	0.1	10
tv-contact-rate-50_4	0.1	10
tv-contact-rate-60_4	0.1	10
tv-contact-rate-70_4	0.1	10
tv-contact-rate-80_4	0.1	10
tv-contact-rate-endday_1	100	190
tv-contact-rate-endday_2	230	305
tv-contact-rate-endday_3	330	420
 	 	 
Constant ODESIM Params	value	 
symp-frac-davies		 
steps-per-day	2	 
time-symp-to-hosp	4	 
tv-vaccinees-endday	364 371 378 385 392 399 406 413 420 427 434 441 448 455 462 469 476 483 490 497 504 511 518 525 532	 
tv-vaccinees-10	0 1 73 50 260 316 594 761 1245 933 1067 1000 1400 1300 2200 5500 7000 12900 15400 23000 25200 25000 19000 65400 56800 0	 
tv-vaccinees-20	0 114 5089 3414 5890 8104 11260 13313 14618 12837 11361 9700 7100 8800 17600 34900 32800 45300 47600 64200 51800 51900 38000 34200 20000 0	 
tv-vaccinees-30	0 262.30 7295.87 4262.15 7001.69 10440.36 12856.47 16457.12 17931.13 14963.25 14600.10 13794.38 11861.08 16354.71 29626.57 49377.62 39136.34 49273.12 49847.89 65732.33 47966.83 47548.82 33963.45 30096.84 16615.96 0	 
tv-vaccinees-40	0 239.70 6667.13 3894.85 6398.31 9540.64 11748.53 15038.88 16385.87 13673.75 13341.90 12605.62 10838.92 14945.29 27073.43 45122.38 35763.66 45026.88 45552.11 60067.67 43833.17 43451.18 31036.55 27503.16 15184.04 0	 
tv-vaccinees-50	0 302.36 6421.90 3455.52 7674.71 11654.15 16508.68 17385.39 20844.83 16785.03 20586.38 22143.13 22950.78 31633.04 45026.61 67573.56 65285.21 74640.52 65419.82 69727.30 39171.13 33517.56 25710.26 21537.39 11980.17 0	 
tv-vaccinees-60	0 262.64 5578.10 3001.48 6666.29 10122.85 7400.79 13304.32 17623.59 15180.12 27705.30 36443.52 48511.61 66421.82 64171.44 61272.81 52890.35 51597.17 40961.56 41090.43 23805.31 20793.52 18115.33 13965.33 7942.69 0	 
tv-vaccinees-70	0 18.13 377.76 283.61 3830.21 4627.25 7712.89 7061.97 13743.49 29800.61 45812.12 45011.99 42377.18 49766.45 41376.92 30098.11 21665.27 16083.03 9896.04 8095.06 5516.23 5216.96 6073.19 4099.10 2491.16 0	 
tv-vaccinees-80	0 12.87 268.24 201.39 2719.79 3285.75 9027.64 5792.32 14899.09 49401.24 63159.20 46701.36 19560.43 9678.69 8325.03 10355.52 6159.17 5279.28 3722.57 3587.21 2707.33 2571.96 2301.23 2098.18 1285.98 0	 
introday	55	 
