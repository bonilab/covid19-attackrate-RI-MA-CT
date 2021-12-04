#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
//#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "essentials.h"
#include "prms.h"
#include "assert.h"


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern FILE* OutFile;
extern double* yic;
extern double* s_to_r;
extern prms* ppc;

extern double G_C_COMM[NUMAC][NUMAC];

extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;

// extern int G_CLO_FIRSTLOCKDOWN_ENDDAY;

extern bool G_B_DIAGNOSTIC_MODE;
extern bool G_B_CHECKPOP_MODE;
extern bool G_B_REMOVE_99COLUMNS;
extern bool G_B_BINARY_OUTPUT;

extern bool G_B_S_TO_R_OUTPUT;

// extern double G_CLO_ICUFRAC_DEV;
// extern double G_CLO_ICUFRAC_DEV_SECONDPHASE;
// extern int G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY;

extern int G_CLO_STEPS_PER_DAY;


//  END  ### ### GLOBAL VARIABLES ### ###

void write(double tt, double * yic, size_t dim)
{
    if( G_B_BINARY_OUTPUT){
        write_bin(tt, yic, dim);
        return;
    }
    int startcol = 0;
    if( G_B_REMOVE_99COLUMNS ){ startcol=99; }

    if( OutFile == NULL )
    {
        printf("%1.3f", tt);
        for(int i=startcol; i<dim; i++) printf("\t%1.4f", yic[i]);
        //for(i=0;i<NUMAC+3;i++) printf("\t%1.4f", yic[i]);
        printf("\n");
    }
    else
    {
        fprintf(OutFile, "%1.3f", tt);
        for(int i=startcol; i<dim; i++) fprintf(OutFile, "\t%1.4f", yic[i]);
        //for(i=0;i<NUMAC+3;i++) printf("\t%1.4f", yic[i]);
        fprintf(OutFile, "\n");
    }
}

void write_bin(double tt, double * yic, size_t dim)
{
    int startcol = 0;
    if( G_B_REMOVE_99COLUMNS ){ startcol=99; }

    if( OutFile == NULL )
    {
        fwrite(&tt, sizeof(tt), 1, stdout);
        fwrite(yic+startcol, sizeof(double), dim-startcol, stdout);
    }
    else
    {
        fwrite(&tt, sizeof(tt), 1, OutFile);
        fwrite(yic+startcol, sizeof(double), dim-startcol, OutFile);
    }
}

// [[Rcpp::export]]
void generate_trajectories( double inital_number_of_cases, double param_beta, double t0, double tf, double h )
{
    assert(ppc);
    double tt,rkstep,ttstop,eps,h1,hmin,ttbeforestop;
    int nvar,nok,nbad; //rkqs();
    int i, j;
    int dim = STARTK+NUMAC; // this is the dimension of the dynamical system, i.e. the number of equations (May 20 2020: this should be 306)
    
    int startcol=0;         // this is the starting column of the output; since the first 99 columns (indices 0 to 98) are the S, E, A classes,
                            // these may not need to be output since they don't have likelihoods associated with them
                            
    if( G_B_REMOVE_99COLUMNS ) startcol=99;
    
    double totalpop=0.0;
    if( G_B_CHECKPOP_MODE )
    {
        for(int k=0; k<dim-18; k++) totalpop += yic[k];
        printf("\n\t totalpop=%1.16f", totalpop);
    }
    
    
    // NOTE 2020-04-04 : this is still very fast and it prevents an off-by-one-day pseudo-error when the step size is close to 0.5 or 1.0
    //                 : just to be clear, a step-size of 1.0 is perfectly fine, but you have to mentally correct for the fact that the 
    //                 : diff-eqs will be about one day late
    int steps_per_day = G_CLO_STEPS_PER_DAY;
    
    // set some values in the RKF integrator
    tt = (double)((int)t0); // chop off any decimals, and just start at an integer time
    ttstop = (double)((int)tf); // chop off any decimals, and just start at an integer time
    rkstep = 1.0 / ((double)steps_per_day);
    eps=0.00000001;
    h1 = 0.1;
    hmin = 1e-13;
    //hmin =    0.00000000001;
    nvar = dim;
    
    int counter = 0;
    bool bIntroductionOccurred = false;
    // bool bContactRatesUpdatedAfterLockdown = false;
    bool bPhase2BetterClinicalManagementBegun = false;
    double NextTimeForBetaUpdate=1000000.0; assert( tf < 999999.0 );
    
    // if the first tv-contact-rate-beginday CLO is higher than 0, use the default contact rates (1.0 for all groups) first 
    bool using_default_contact_rate_first = (ppc->tvp_contact_rate->get_day_index() == 0) && (ppc->tvp_contact_rate->get_current_begin_day() > 0);
    // update at the end of end_day
    double next_end_day_tvp_contact_rate = using_default_contact_rate_first ? ((double)ppc->tvp_contact_rate->get_current_begin_day()) - 0.000001 : ((double)ppc->tvp_contact_rate->get_current_end_day()) + 0.999999; 
    bool end_tvp_contact_rate = false;

    // if the first tv-contact-coeff-beginday CLO is higher than 0, use the default contact rates (1.0 for all groups) first 
    bool using_default_contact_coeff_first = (ppc->tvp_contact_coeff->get_day_index() == 0) && (ppc->tvp_contact_coeff->get_current_begin_day() > 0);
    // update at the end of end_day
    double next_end_day_tvp_contact_coeff = using_default_contact_coeff_first ? ((double)ppc->tvp_contact_coeff->get_current_begin_day()) - 0.000001 : ((double)ppc->tvp_contact_coeff->get_current_end_day()) + 0.999999; 
    bool end_tvp_contact_coeff = false;

    bool using_contac_rate = !ppc->tvp_contact_rate->empty_data() && ppc->tvp_contact_coeff->empty_data(); 

    // if the first tv-contact-rate-beginday CLO is higher than 0, use the default contact rates (1.0 for all groups) first 
    bool using_default_hosp_frac_first = (ppc->tvp_hosp_frac->get_day_index() == 0) && (ppc->tvp_hosp_frac->get_current_begin_day() > 0);
    // update at the end of end_day
    double next_end_day_tvp_hosp_frac = using_default_hosp_frac_first ? ((double)ppc->tvp_hosp_frac->get_current_begin_day()) - 0.000001 : ((double)ppc->tvp_hosp_frac->get_current_end_day()) + 0.999999; 
    bool end_tvp_hosp_frac = false;

    // vaccination
    double next_start_day_tvp_vac_num = ((double)ppc->tvp_vac_num->get_current_begin_day()) - rkstep + 0.000001 ; 
    double next_end_day_tvp_vac_num = ((double)ppc->tvp_vac_num->get_current_end_day()); 
    double next_num_day_tvp_vac_num = ((double)ppc->tvp_vac_num->get_current_end_day()) - ((double)ppc->tvp_vac_num->get_current_begin_day()) + 1;
    // bool done_vac = ppc->tvp_vac_num->get_current_begin_day() + 1 != ppc->tvp_vac_num->get_current_end_day();
    bool end_tvp_vac_num = ppc->tvp_vac_num->empty_data();

    // if the first tv-dev-icu-frac-beginday CLO is higher than 0, use the default dev-icu-frac of 1.0 first 
    bool using_default_dev_icu_frac_first = (ppc->tvp_dev_icu_frac->get_day_index() == 0) && (ppc->tvp_dev_icu_frac->get_current_begin_day() > 0);
    // update at the end of end_day
    double next_end_day_tvp_dev_icu_frac = using_default_dev_icu_frac_first ? ((double)ppc->tvp_dev_icu_frac->get_current_begin_day()) - 0.000001 : ((double)ppc->tvp_dev_icu_frac->get_current_end_day()) + 0.999999; 
    bool end_tvp_dev_icu_frac = false;

    // if the first tv-dev-len-hospstay-beginday CLO is higher than 0, use the default dev-len-hospstay of 1.0 first 
    bool using_default_dev_len_hospstay_first = (ppc->tvp_dev_len_hospstay->get_day_index() == 0) && (ppc->tvp_dev_len_hospstay->get_current_begin_day() > 0);
    // update at the end of end_day
    double next_end_day_tvp_dev_len_hospstay = using_default_dev_len_hospstay_first ? ((double)ppc->tvp_dev_len_hospstay->get_current_begin_day()) - 0.000001 : ((double)ppc->tvp_dev_len_hospstay->get_current_end_day()) + 0.999999; 
    bool end_tvp_dev_len_hospstay = false;
    

    if( ppc->v_betatimes.size() > 1 )
    {
        NextTimeForBetaUpdate = ppc->v_betatimes[1];
    }
    ppc->assign_new_beta(); // called here for the first time, so it just puts the first beta into use
    
    // before integration begins, assign the initial beta value
    
    //
    //BEGIN OF LOOP INTEGRATING ODEs
    //
    while( tt < (ttstop-0.000000001) )      // MFB NOTE 2020-06-09: subtract this tiny amount since sometimes we shouldn't enter the while loop
    {                                       // but we do anyway because tt is just barely below ttstop by some tiny fraction
        // introduce the first infections
        if( !bIntroductionOccurred && tt > G_CLO_INTRODUCTION_TIME )
        {
            yic[STARTI+4]=((double)G_CLO_INTRODUCTION_COUNT); // some number of infected individuals in their forties are introduced 
            bIntroductionOccurred = true;
        }
        
        // if( !bContactRatesUpdatedAfterLockdown && tt > G_CLO_FIRSTLOCKDOWN_ENDDAY )
        if ( using_contac_rate ){
            if( tt > next_end_day_tvp_contact_rate && tt < next_end_day_tvp_contact_rate + rkstep && !end_tvp_contact_rate )
            {
                // if not using default contact rate first (because the first tv-contact-rate-beginday CLO is 0), update indices
                // else, keep the current indices (which should be 0)
                if (!using_default_contact_rate_first){
                    end_tvp_contact_rate = !ppc->tvp_contact_rate->advance_indices(); 
                } else {
                    using_default_contact_rate_first = false; // reset so that advance_indices() can be called next time
                }
                            
                // printf("\n#### contact-rate state ####\ntime %.3f\n", tt);
                // ppc->tvp_contact_rate->print_current_state();

                for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
                {
                    for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
                    {
                        // G_C_COMM[acs][aci] = ppc->v_mixing_level_postld[acs] * ppc->v_mixing_level_postld[aci];
                        G_C_COMM[acs][aci] = ppc->tvp_contact_rate->get_current_data(acs) * ppc->tvp_contact_rate->get_current_data(aci);
                    }
                }        
                // bContactRatesUpdatedAfterLockdown = true;
                next_end_day_tvp_contact_rate = (double)ppc->tvp_contact_rate->get_current_end_day() + 0.999999; // update at the end of end_day
            }
        } else
        {
            if( tt > next_end_day_tvp_contact_coeff && tt < next_end_day_tvp_contact_coeff + rkstep && !end_tvp_contact_coeff )
            {
                // if not using default contact rate first (because the first tv-contact-rate-beginday CLO is 0), update indices
                // else, keep the current indices (which should be 0)
                if (!using_default_contact_coeff_first){
                    end_tvp_contact_coeff = !ppc->tvp_contact_coeff->advance_indices(); 
                } else {
                    using_default_contact_coeff_first = false; // reset so that advance_indices() can be called next time
                }
                            
                // printf("\n#### contact-rate state ####\ntime %.3f\n", tt);
                // ppc->tvp_contact_coeff->print_current_state();

                if( ppc->tvp_contact_coeff->get_data_index() >= 1){
                    // CoMix - Belgium wave 5 (late June 2020) matrix from `socialmixr` R package;
                    G_C_COMM[0][0] = 0.2;
                    G_C_COMM[0][1] = 0.136363636363636;     G_C_COMM[1][0] = 0.136363636363636; 
                    G_C_COMM[0][2] = 0.165178571428571;     G_C_COMM[2][0] = 0.165178571428571;
                    G_C_COMM[0][3] = 1.11034482758621;      G_C_COMM[3][0] = 1.11034482758621;
                    G_C_COMM[0][4] = 0.313218390804598;     G_C_COMM[4][0] = 0.313218390804598;
                    G_C_COMM[0][5] = 0.072769953051643;     G_C_COMM[5][0] = 0.072769953051643;
                    G_C_COMM[0][6] = 0.135514018691589;     G_C_COMM[6][0] = 0.135514018691589;
                    G_C_COMM[0][7] = 0.048076923076923;     G_C_COMM[7][0] = 0.048076923076923;
                    G_C_COMM[0][8] = 0.001;                 G_C_COMM[8][0] = 0.001;

                    G_C_COMM[1][1] = 1.13636363636364;
                    G_C_COMM[1][2] = 0.001;                 G_C_COMM[2][1] = 0.59375;
                    G_C_COMM[1][3] = 1.13636363636364;      G_C_COMM[3][1] = 0.293103448275862;
                    G_C_COMM[1][4] = 0.863636363636364;     G_C_COMM[4][1] = 0.844827586206897;
                    G_C_COMM[1][5] = 0.001;                 G_C_COMM[5][1] = 0.375586854460094;
                    G_C_COMM[1][6] = 0.001;                 G_C_COMM[6][1] = 0.226635514018692;
                    G_C_COMM[1][7] = 0.001;                 G_C_COMM[7][1] = 0.173076923076923;
                    G_C_COMM[1][8] = 0.001;                 G_C_COMM[8][1] = 0.001;

                    G_C_COMM[2][2] = 2.61607142857143;
                    G_C_COMM[2][3] = 1.08035714285714;      G_C_COMM[3][2] = 1.16896551724138;
                    G_C_COMM[2][4] = 1.29464285714286;      G_C_COMM[4][2] = 0.985632183908046;
                    G_C_COMM[2][5] = 0.763392857142857;     G_C_COMM[5][2] = 1.22300469483568;
                    G_C_COMM[2][6] = 0.044642857142857;     G_C_COMM[6][2] = 0.628504672897196;
                    G_C_COMM[2][7] = 0.0625;                G_C_COMM[7][2] = 0.336538461538462;
                    G_C_COMM[2][8] = 0.15625;               G_C_COMM[8][2] = 0.001;

                    G_C_COMM[3][3] = 1.3448275862069;
                    G_C_COMM[3][4] = 0.389655172413793;     G_C_COMM[4][3] = 1.12931034482759;
                    G_C_COMM[3][5] = 0.513793103448276;     G_C_COMM[5][3] = 0.434272300469484;
                    G_C_COMM[3][6] = 0.16551724137931;      G_C_COMM[6][3] = 0.836448598130841;
                    G_C_COMM[3][7] = 0.037931034482759;     G_C_COMM[7][3] = 0.403846153846154;
                    G_C_COMM[3][8] = 0.010344827586207;     G_C_COMM[8][3] = 0.001;

                    G_C_COMM[4][4] = 1.02298850574713;
                    G_C_COMM[4][5] = 0.195402298850575;     G_C_COMM[5][4] = 0.910798122065728;
                    G_C_COMM[4][6] = 0.086206896551724;     G_C_COMM[6][4] = 0.380841121495327;
                    G_C_COMM[4][7] = 0.21551724137931;      G_C_COMM[7][4] = 0.322115384615385;
                    G_C_COMM[4][8] = 0.0001;                G_C_COMM[8][4] = 0.6;

                    G_C_COMM[5][5] = 0.711267605633803;
                    G_C_COMM[5][6] = 0.119718309859155;     G_C_COMM[6][5] = 0.647196261682243;
                    G_C_COMM[5][7] = 0.150234741784038;     G_C_COMM[7][5] = 0.173076923076923;
                    G_C_COMM[5][8] = 0.096244131455399;     G_C_COMM[8][5] = 0.001;

                    G_C_COMM[6][6] = 0.457943925233645;
                    G_C_COMM[6][7] = 0.294392523364486;     G_C_COMM[7][6] = 0.605769230769231;
                    G_C_COMM[6][8] = 0.135514018691589;     G_C_COMM[8][6] = 0.6;

                    G_C_COMM[7][7] = 1.27884615384615;
                    G_C_COMM[7][8] = 0.197115384615385;     G_C_COMM[8][7] = 0.2;

                    G_C_COMM[8][8] = 0.8;  
                }

                for(int aci=0; aci<NUMAC; aci++){
                    G_C_COMM[0][aci] *= ppc->tvp_contact_coeff->get_current_data(0);
                    G_C_COMM[1][aci] *= ppc->tvp_contact_coeff->get_current_data(1);
                    G_C_COMM[2][aci] *= ppc->tvp_contact_coeff->get_current_data(2);
                    G_C_COMM[3][aci] *= ppc->tvp_contact_coeff->get_current_data(3);
                    G_C_COMM[4][aci] *= ppc->tvp_contact_coeff->get_current_data(4);
                    G_C_COMM[5][aci] *= ppc->tvp_contact_coeff->get_current_data(5);
                    G_C_COMM[6][aci] *= ppc->tvp_contact_coeff->get_current_data(6);
                    G_C_COMM[7][aci] *= ppc->tvp_contact_coeff->get_current_data(7);
                    G_C_COMM[8][aci] *= ppc->tvp_contact_coeff->get_current_data(8);
                }         
                // bContactRatesUpdatedAfterLockdown = true;
                next_end_day_tvp_contact_coeff = (double)ppc->tvp_contact_coeff->get_current_end_day() + 0.999999; // update at the end of end_day
            }
        }
        

        // time-varying hosp frac
        if( tt > next_end_day_tvp_hosp_frac && tt < next_end_day_tvp_hosp_frac + rkstep && !end_tvp_hosp_frac )
        {
            // if not using default contact rate first (because the first tv-contact-rate-beginday CLO is 0), update indices
            // else, keep the current indices (which should be 0)
            if (!using_default_hosp_frac_first){
                for (int ac = 0; ac < NUMAC; ac++){
                    ppc->v_prob_I2_H[ac] = ppc->tvp_hosp_frac->get_current_data(ac);
                    // fprintf(stdout, "v_prob_I2_H[%d] %f\n", ac, ppc->v_prob_I2_H[ac]);
                    // printf("v_prob_I2_H[%d] %f\n", ac, ppc->v_prob_I2_H[ac]);
                }                
                end_tvp_hosp_frac = !ppc->tvp_hosp_frac->advance_indices(); 
            } else {
                ppc->v_prob_I2_H[0] = 0.021;    // NOTE because there is not likely to be data on hosp rates for 0-9 y.o., we adopt the hosp rate for 10-19 year olds here
                ppc->v_prob_I2_H[1] = 0.021;
                ppc->v_prob_I2_H[2] = 0.017;
                ppc->v_prob_I2_H[3] = 0.031;    // this is the hospitalization probability for 30-39 year-olds
                ppc->v_prob_I2_H[4] = 0.05;
                ppc->v_prob_I2_H[5] = 0.1;
                ppc->v_prob_I2_H[6] = 0.15;
                ppc->v_prob_I2_H[7] = 0.15;
                ppc->v_prob_I2_H[8] = 0.15;
                using_default_hosp_frac_first = false; // reset so that advance_indices() can be called next time
            }
            // fprintf(stdout, "v_prob_I2_H[8] %f\n", ppc->v_prob_I2_H[8]);
                        
            // printf("\n#### contact-rate state ####\ntime %.3f\n", tt);
            // ppc->tvp_hosp_frac->print_current_state();
            next_end_day_tvp_hosp_frac = (double)ppc->tvp_hosp_frac->get_current_end_day() + 0.999999; // update at the end of end_day
        }


        // vaccination
        // if (tt > next_start_day_tvp_vac_num && tt < next_start_day_tvp_vac_num + rkstep && !end_tvp_vac_num && !done_vac)
        if (tt > next_start_day_tvp_vac_num && tt < next_end_day_tvp_vac_num  + 0.999999 && !end_tvp_vac_num ){
            double num_vaccinees;
            
            // fprintf(stdout, "\ntime %f\n", tt);

            for (int ac = 0; ac < NUMAC; ac++){
                num_vaccinees = floor(ppc->tvp_vac_num->get_current_data(ac) / (next_num_day_tvp_vac_num * steps_per_day));
                if (tt > next_end_day_tvp_vac_num + rkstep - 0.000001){
                    num_vaccinees +=  ppc->tvp_vac_num->get_current_data(ac) - num_vaccinees * next_num_day_tvp_vac_num * steps_per_day;
                    // fprintf(stdout, "%f of %f \n", num_vaccinees, ppc->tvp_vac_num->get_current_data(ac));
                }
                ppc->v_cumul_vac[ac] += num_vaccinees; // cumulative vaccinations

                int i = 0;
                // S / (S + E + A + (1 - (1 - prob_EA)*reporting_rate - prob_vaccination)*R)
                //
                // Nov 2021 modification:
                // S / {S + E + A + [1 - (1 - prob_EA)*reporting_rate]*(1 - prob_vaccination)*R }
                //
                double frac_vac_seroneg = yic[ac]; // S
                for (i = 0; i < NUME; i++){ frac_vac_seroneg += yic[STARTE + i*NUMAC + ac]; } // S + E
                for (i = 0; i < NUMA; i++){ frac_vac_seroneg += yic[STARTA + i*NUMAC + ac]; } // S + E + A

                // tmp = = 1.0 - (1.0 - ppc->v_prob_E_A[ac])*ppc->v[i_avg_reporting_rate]*(1.0 - ppc->v_cumul_vac[ac]/ppc->v_pop_ac[ac]); // 
                double tmp = (1.0 - ppc->v_prob_E_A[ac]); // (1 - prob_EA)
                tmp *= ppc->v[i_avg_reporting_rate]; // (1 - prob_EA)*reporting_rate
                tmp = 1.0 - tmp; // 1 - (1 - prob_EA)*reporting_rate
                tmp *= (1.0 - ppc->v_cumul_vac[ac]/ppc->v_pop_ac[ac]); // [1 - (1 - prob_EA)*reporting_rate] * (1 - prob_vaccination)
                tmp = (tmp < 0.0) ? 0.0 : tmp;
                tmp *= yic[STARTR + ac]; // [1 - (1 - prob_EA)*reporting_rate] * (1 - prob_vaccination) * R
                
                frac_vac_seroneg += tmp; // S + E + A + [1 - (1 - prob_EA)*reporting_rate] * (1 - prob_vaccination) * R
                frac_vac_seroneg = yic[ac] / frac_vac_seroneg; // S / {S + E + A + [1 - (1 - prob_EA)*reporting_rate] * (1 - prob_vaccination) * R}

                num_vaccinees *= frac_vac_seroneg;
                num_vaccinees = (num_vaccinees < yic[ac]) ? num_vaccinees : yic[ac];
                
                s_to_r[ac] += num_vaccinees;
                // s_to_r[ac] = frac_vac_seroneg;

                yic[ac] -= num_vaccinees;
                yic[STARTR + ac] += num_vaccinees;

                // fprintf(stdout, "vaccinate %d %f\n", ac, num_vaccinees);
            }
            
            if (tt > next_end_day_tvp_vac_num + rkstep - 0.000001){
                end_tvp_vac_num = !ppc->tvp_vac_num->advance_indices();
                next_start_day_tvp_vac_num = ((double)ppc->tvp_hosp_frac->get_current_begin_day()) - rkstep + 0.000001 ; 
                next_end_day_tvp_vac_num = ((double)ppc->tvp_vac_num->get_current_end_day()); 
                next_num_day_tvp_vac_num = ((double)ppc->tvp_vac_num->get_current_end_day()) - ((double)ppc->tvp_vac_num->get_current_begin_day()) + 1;
                // done_vac = ppc->tvp_vac_num->get_current_begin_day() + 1 != ppc->tvp_vac_num->get_current_end_day();
            }
        }
        


        /* if( !bPhase2BetterClinicalManagementBegun && tt >= G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY )
        {
            for(int ac=0; ac<NUMAC; ac++)
            {
                ppc->v_prob_HA_CA[ac] *= (G_CLO_ICUFRAC_DEV_SECONDPHASE / G_CLO_ICUFRAC_DEV);
            }            

            bPhase2BetterClinicalManagementBegun = true;
        } */
        if( tt > next_end_day_tvp_dev_icu_frac && tt < next_end_day_tvp_dev_icu_frac + rkstep && !end_tvp_dev_icu_frac )
        {
            double dif_prev, dif_curr;
            // if not using default contact rate first (because the first tv-dev-icu-frac-beginday CLO is 0), update indices
            // else, keep the current indices (which should be 0)
            if (!using_default_dev_icu_frac_first){
                dif_prev = ppc->tvp_dev_icu_frac->get_current_data();
                end_tvp_dev_icu_frac = !ppc->tvp_dev_icu_frac->advance_indices(); 
            } else {
                dif_prev = 1.0;
                using_default_dev_icu_frac_first = false; // reset so that advance_indices() can be called next time
            }
            dif_curr = ppc->tvp_dev_icu_frac->get_current_data() / dif_prev;
                        
            // printf("\n#### dev-icu-frac state ####\ntime %.3f\n", tt);
            // ppc->tvp_dev_icu_frac->print_current_state();

            for(int ac=0; ac<NUMAC; ac++) // age-class of the susceptible individual
            {
                // ppc->v_prob_HA_CA[ac] *= (G_CLO_ICUFRAC_DEV_SECONDPHASE / G_CLO_ICUFRAC_DEV);
                ppc->v_prob_HA_CA[ac] *= dif_curr;// / dif_prev;
            }        
            // bPhase2BetterClinicalManagementBegun = true;
            next_end_day_tvp_dev_icu_frac = (double)ppc->tvp_dev_icu_frac->get_current_end_day() + 0.999999; // update at the end of end_day
        }
        
        
        if( tt > next_end_day_tvp_dev_len_hospstay && tt < next_end_day_tvp_dev_len_hospstay + rkstep && !end_tvp_dev_len_hospstay )
        {
            double dlh_prev, dlh_curr;
            // if not using default contact rate first (because the first tv-dev-len-hospstay-beginday CLO is 0), update indices
            // else, keep the current indices (which should be 0)
            // ppc->v[ i_len_medicalfloor_hospital_stay ]              = 10.7 * G_CLO_DEV_LEN_HOSPSTAY;
            // fprintf(stdout, "\nppc->v[ i_len_medicalfloor_hospital_stay ] %f\n", ppc->v[ i_len_medicalfloor_hospital_stay ]);
            if (!using_default_dev_len_hospstay_first){
                dlh_prev = ppc->tvp_dev_len_hospstay->get_current_data();
                // ppc->v[ i_len_medicalfloor_hospital_stay ] = 10.7 * ppc->tvp_dev_len_hospstay->get_current_data();
                end_tvp_dev_len_hospstay = !ppc->tvp_dev_len_hospstay->advance_indices(); 
            } else {
                dlh_prev = 1.0;
                // ppc->v[ i_len_medicalfloor_hospital_stay ] = 10.7 ;
                using_default_dev_len_hospstay_first = false; // reset so that advance_indices() can be called next time
            }
            dlh_curr = ppc->tvp_dev_len_hospstay->get_current_data() / dlh_prev;
            ppc->v[ i_len_medicalfloor_hospital_stay ] *= dlh_curr;

            // printf("\n#### tvp_dev_len_hospstay state ####\ntime %.3f\n", tt);
            // ppc->tvp_dev_len_hospstay->print_current_state();
            // fprintf(stdout, "\nppc->v[ i_len_medicalfloor_hospital_stay ] %f\n", ppc->v[ i_len_medicalfloor_hospital_stay ]);

            next_end_day_tvp_dev_len_hospstay = (double)ppc->tvp_dev_len_hospstay->get_current_end_day() + 0.999999; // update at the end of end_day
        }
        
        
        // check if the hospitalization rates need to be updated (because we are out of the early period)
        /*if( ppc->earlymarch_highhosp_period )
        {
            if( tt > ppc->earlymarch_highhosp_endday )
            {
                ppc->end_earlymarch_hosprates();
            }
        }*/
        
        // check if the beta value needs to be updated
        if( tt > NextTimeForBetaUpdate )
        {
            ppc->assign_new_beta();
            NextTimeForBetaUpdate = ppc->get_new_update_time();
        }
        


        //
        ////////////// integrate
        //
        if ( odeint(yic,nvar,tt,tt+rkstep,eps,&h1,hmin,&nok,&nbad, derivs,rkqs) != 1)
      	{
	        fprintf(stderr, "\n\nEXITING BC OF ERROR IN ODEINT\n\n");
            exit(-1);
        }




        
        if( counter%steps_per_day==0 && !G_B_DIAGNOSTIC_MODE && !G_B_CHECKPOP_MODE && tt>49.5 )
        {
            if ( G_B_S_TO_R_OUTPUT ){
                write(tt, s_to_r, NUMAC);
                // write(tt, &ppc->v_cumul_vac[0], NUMAC);
            } else{
                write(tt, yic, dim);
            }
        }
        
        // ### BEGIN only enter block below if you're checking the total population size
        if( counter%(50*steps_per_day)==0 && G_B_CHECKPOP_MODE )
        {
            totalpop=0.0;
            for(int k=0; k<dim-18; k++) totalpop += yic[k];
            printf("\n\t totalpop=%1.16f", totalpop);
            
        }
        // ### END only enter block if you're checking the total population size
        
        tt += rkstep;
        counter++;
    }
    //
    //END OF LOOP INTEGRATING ODEs
    //

    // print out results of the last step
    if( !G_B_DIAGNOSTIC_MODE && !G_B_CHECKPOP_MODE )
    {
        if ( G_B_S_TO_R_OUTPUT ){
            write(tt, s_to_r, NUMAC);
            // write(tt, &ppc->v_cumul_vac[0], NUMAC);
        } else{
            write(tt, yic, dim);
        }
    }

    if( G_B_CHECKPOP_MODE ) printf("\n\n");
    //printf("\n\n    %d \n\n", Q );
    
    //delete[] yic;    
    
    //return vv;
}

