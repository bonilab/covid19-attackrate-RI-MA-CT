#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "prms.h"
#include "essentials.h"
#include "parseargs.h"

#include <ctime>
#include <sys/time.h>

using namespace std;


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double b,d,s,c0,c2,c3;	// this are defined, i believe, in rkf.cpp
double* yic;
double* s_to_r; // vaccinees who are in S

FILE* OutFile = NULL;

prms* ppc;  // need to allocate space for this in main

double G_CLO_BETA = 0.3;
double G_CLO_INTRODUCTION_TIME = -1;
int G_CLO_INTRODUCTION_COUNT;
string G_CLO_LOCATION = "RI";

double G_CLO_TF = 365.0;
int G_CLO_STEPS_PER_DAY=100;
double G_CLO_P_HOSP_TO_ICU = 0.30;

double G_CLO_SYMP_FRAC = 0.25;                //SHOULD BE DEPRECATED SOON

double G_CLO_SYMP_FRAC_EQUAL = 0.601;           // this number should be = 0.564; it is taken from Fuhan's lit review July 8 2020;
                                                // if we exclude the studies that did not distinguish btw pre-symptomatic and asymptomatic, this number is = 0.601

double G_B_SYMP_FRAC_DAVIES_20200616 = false;
double G_B_SYMP_FRAC_SIMPLEAVERAGE   = true;           // this comes from simple weighted averages within 20-yr age groups, from Fuhan's lit review July 8 2020


//double G_CLO_HOSPFRAC_YOUNG_DEV = 1.0;      //DEPRECATED
//double G_CLO_HOSPFRAC_MID_DEV = 0.5;        //DEPRECATED
//double G_CLO_HOSPFRAC_OLD_DEV = 0.5;        //DEPRECATED

/* double G_CLO_HOSPFRAC_10 = 0.021; // numbers below are the raw RI data (June 8)
double G_CLO_HOSPFRAC_20 = 0.017;
double G_CLO_HOSPFRAC_30 = 0.031;
double G_CLO_HOSPFRAC_40 = 0.05;
double G_CLO_HOSPFRAC_50 = 0.1;          //0.125;
double G_CLO_HOSPFRAC_60 = 0.15;         //0.222;
double G_CLO_HOSPFRAC_70 = 0.15;         //0.307;
double G_CLO_HOSPFRAC_80 = 0.15;         //0.198; */

// double G_CLO_MIXINGLEVEL[NUMAC];
// double G_CLO_MIXINGLEVEL_POSTLD[NUMAC];     // the age-specific mixing level post-lockdown
// int G_CLO_FIRSTLOCKDOWN_ENDDAY = 121;       // this is April 30 2020
// bool G_CLO_POSTLD_MIXING_SET = false;

double G_CLO_DEATHPROB_HOME_60 = 0.00;
double G_CLO_DEATHPROB_HOME_70 = 0.05;
double G_CLO_DEATHPROB_HOME_80 = 0.22;

double G_CLO_DEATHPROB_POSTVENT = 0.189; // this is the average death probablity for ICU but extubated patients, which applies roughly to patients 60 and over (no data for younger patients)


double G_CLO_PROB_NONICU_DEATH_80 = 0.05; // probability of dying in the hospital but not in the critical care unit

double G_CLO_MIN_MIXINGLEVEL_00 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_10 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_20 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_30 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_40 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_50 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_60 = 0.00;
double G_CLO_MIN_MIXINGLEVEL_70 = 0.00;     // the mininum ("floor") relative mixing level that certain age groups can 
double G_CLO_MIN_MIXINGLEVEL_80 = 0.00;     // go down to; if you set this to 0.15, then people in that age group will be able to
                                            // reduce their mixing level by 85% during a lockdown, but no more than that


double G_CLO_RELSUSC_0 = 0.6;
double G_CLO_RELSUSC_10 = 0.6;
double G_CLO_RELSUSC_20 = 1.0;
double G_CLO_RELSUSC_30 = 1.0;
double G_CLO_RELSUSC_40 = 1.0;
double G_CLO_RELSUSC_50 = 1.0;
double G_CLO_RELSUSC_60 = 1.0;
double G_CLO_RELSUSC_70 = 1.0;
double G_CLO_RELSUSC_80 = 1.0;

double G_CLO_TIME_SYMPTOHOSP = 7.0;
double G_CLO_SELFISOLATION_FACTOR = 1.0;

double G_CLO_EARLYMARCH_HOSPRATE = 1.0;     // scaling factor showing how much more likely hospitalization was in early March than later in the epidemic
double G_CLO_EARLYMARCH_ENDDAY = 84.5;      // the day that the high-hosp rate ended (probably because the patient load became too high)


// double G_CLO_DEV_LEN_HOSPSTAY = 1.0;

// double G_CLO_ICUFRAC_DEV = 1.0;
// double G_CLO_ICUFRAC_DEV_SECONDPHASE = 1.0;             // the second phase of the epidemic, when clinical management has improved, and
//                                                         // and fewer patients move to the ICU
// int G_CLO_BETTERCLINICALMANAGEMENT_BEGINDAY = 153;      // this is June 1 2020

double G_CLO_PROB_ICU_TO_VENT = 0.75;

double G_CLO_VENTDEATH_MID_DEV = 0.7;
double G_CLO_VENTDEATH_70_DEV = 1.0;
double G_CLO_VENTDEATH_80_DEV = 1.0;

double G_CLO_MEANTIME_ON_VENT_SURV = 10.8;   //  mean time on a ventilator for survivors

double G_CLO_RELATIVE_BETA_HOSP = 0.2;

bool G_B_DIAGNOSTIC_MODE = false;
bool G_B_CHECKPOP_MODE = false;
bool G_B_USESOCIALCONTACTMATRIX = true;
bool G_B_REMOVE_99COLUMNS = false;
bool G_B_BINARY_OUTPUT = false;
bool G_B_DEATHRATE_OUTPUT = false;

bool G_B_S_TO_R_OUTPUT = false;

double G_C_COMM[NUMAC][NUMAC]; // this is the contact rate in the community between susc (first index) and infected (second index)                              - this matrix is symmetric
double G_C_HOSP[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected hospitalized patient (second index)     - this matrix is not symmetric
double G_C_ICU[NUMAC][NUMAC];  // this is the contact rate btw the community (first index, susceptible) and an infected ICU patient (second index)              - this matrix is not symmetric
double G_C_VENT[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected Ventilated patient (second index)       - this matrix is not symmetric
 

// END ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void InitializeContactMatrices( void );


//
//
int main(int argc, char* argv[])
{
    

    //
    // ###  1.  ALLOCATE SPACE FOR A PARAMETERS CLASS AND FOR THE MAIN STATE VARIABLES
    //
    
    ppc = new prms; assert( ppc );
    yic = new double[STARTK+NUMAC];

    s_to_r = new double[NUMAC];
    for(int i=0; i< NUMAC; i++) s_to_r[i] = 0.0;

    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out

    // assign some default values to the v_beta and v_betatimes arrays
    ppc->v_betas.push_back( 1.0 );
    ppc->v_betas.push_back( 0.9 );
    ppc->v_betas.push_back( 0.8 );
    ppc->v_betatimes.push_back( 0.0 );
    ppc->v_betatimes.push_back( 60.0 );
    ppc->v_betatimes.push_back( 95.0 );
    
    // INITIALIZE GLOBAL VARIABLES
    /* for(int ac=0;ac<NUMAC;ac++) G_CLO_MIXINGLEVEL[ac] = 1.0;
    for(int ac=0;ac<NUMAC;ac++) G_CLO_MIXINGLEVEL_POSTLD[ac] = 1.0; */
    
    
//     string s1(argv[0]); string s2(argv[1]); string s3("3.8.44"); 
//     
//     printf("\n\t %s float outcome %d", s1.c_str() , isFloat(s1)?1:0 );
//     printf("\n\t %s float outcome %d", s2.c_str() , isFloat(s2)?1:0 );
//     printf("\n\t %s float outcome %d", s3.c_str() , isFloat(s3)?1:0 );
//     printf("\n\n");
    
    // this just sets all the matrix values to one
    
    ParseArgs( argc, argv );


    // NOTE you have to assign these mixing levels before you initialize the contact InitializeContactMatrices
    //
    // assign the mixing levels by age group
    /* for(int ac=0;ac<NUMAC;ac++) ppc->v_mixing_level[ac] = G_CLO_MIXINGLEVEL[ac];
    for(int ac=0;ac<NUMAC;ac++) ppc->v_mixing_level_postld[ac] = G_CLO_MIXINGLEVEL_POSTLD[ac];

    if( !G_CLO_POSTLD_MIXING_SET )
    {
        for(int ac=0;ac<NUMAC;ac++) 
        {
            // if the post-LD levels were not set, make sure they are equal to the pre-LD levels
            if( fabs(ppc->v_mixing_level[ac] - ppc->v_mixing_level_postld[ac]) > 0.000001 )
            {
                assert( false );
            }
        }
    } */

    // ppc->tvp_contact_rate->print_all();
    // ppc->tvp_contact_coeff->print_all();

    // if no contact rate is set via command-line options, use default
    if (!ppc->tvp_contact_rate->empty_data() && !ppc->tvp_contact_coeff->empty_data() ) {
        fprintf(stderr,"\n\n\t-tv-contact-rate-* and -tv-contact-coeff-* command-line options are mutually exclusive.\n\n");
        ppc->tvp_contact_rate->print_all();
        ppc->tvp_contact_coeff->print_all();
        exit(-1);
    } 
    if (ppc->tvp_contact_rate->empty_data() && ppc->tvp_contact_coeff->empty_data() ) {
        for (size_t i = 0; i < NUMAC; i++) { ppc->tvp_contact_rate->append_data(i, 1.0); }
    } 
    
    if ( !ppc->tvp_contact_rate->empty_data() ){
        // fill missing begin and/or end days in tvp_contact_rate
        // advance if good_days() and good_integrity()
        // else, print error message
        if ( ppc->tvp_contact_rate->fill_missing_days() ){
            ppc->tvp_contact_rate->advance_indices(); // start day_index and data_index at 0
        } else {
            fprintf(stderr,"\n\n\tError in parsing contact rates.\n\tPlease check all -tv-contact-rate-* command-line options or unset all these arguments to use the default values.\n\n");
            ppc->tvp_contact_rate->print_all();
            exit(-1);
        }
    } else if ( !ppc->tvp_contact_coeff->empty_data() ){
        // fill missing begin and/or end days in tvp_contact_coeff
        // advance if good_days() and good_integrity()
        // else, print error message
        if ( ppc->tvp_contact_coeff->fill_missing_days() ){
            ppc->tvp_contact_coeff->advance_indices(); // start day_index and data_index at 0
        } else {
            fprintf(stderr,"\n\n\tError in parsing contact coefficients.\n\tPlease check all -tv-contact-coeff-* command-line options or unset all these arguments to use the default values.\n\n");
            ppc->tvp_contact_coeff->print_all();
            exit(-1);
        }
    } else {
        fprintf(stderr,"\n\n\tNot enough information to form contact matrices.\n\n");
        exit(-1);
    }
    
    


    // if no hosp frac is set via command-line options, use default (raw RI data as of June 8)
    if (ppc->tvp_hosp_frac->empty_data()) {
        ppc->tvp_hosp_frac->append_data(0, 0.021);
        ppc->tvp_hosp_frac->append_data(1, 0.021);
        ppc->tvp_hosp_frac->append_data(2, 0.017);
        ppc->tvp_hosp_frac->append_data(3, 0.031);
        ppc->tvp_hosp_frac->append_data(4, 0.05);
        ppc->tvp_hosp_frac->append_data(5, 0.1);
        ppc->tvp_hosp_frac->append_data(6, 0.15);
        ppc->tvp_hosp_frac->append_data(7, 0.15);
        ppc->tvp_hosp_frac->append_data(8, 0.15);
    }
    // fill missing begin and/or end days in tvp_hosp_frac
    // advance if good_days() and good_integrity()
    // else, print error message
    if ( ppc->tvp_hosp_frac->fill_missing_days() ){
        ppc->tvp_hosp_frac->advance_indices(); // start day_index and data_index at 0
    } else {
        fprintf(stderr,"\n\n\tError in parsing contact rates.\n\tPlease check all -tv-hosp-frac-* command-line options or unset all these arguments to use the default values.\n\n");
        ppc->tvp_hosp_frac->print_all();
        exit(-1);
    }
    for (int ac = 0; ac < NUMAC; ac++){
        ppc->v_prob_I2_H[ac] = ppc->tvp_hosp_frac->get_current_data(ac);
    } 


    // if no vaccination number is set via command-line options, use default
    if (ppc->tvp_vac_num->empty_data()) {
        // for (size_t i = 0; i < NUMAC; i++) { ppc->tvp_vac_num->append_data(i, 0.0); }
    } else
    {
        // fill missing begin and/or end days in tvp_vac_num
        // advance if good_days() and good_integrity()
        // else, print error message
        if ( ppc->tvp_vac_num->fill_missing_days() ){
            ppc->tvp_vac_num->advance_indices(); // start day_index and data_index at 0
        } else {
            fprintf(stderr,"\n\n\tError in parsing vaccinees.\n\tPlease check all -tv-vaccinees-* command-line options or unset all these arguments to use the default values.\n\n");
            ppc->tvp_vac_num->print_all();
            exit(-1);
        }
    }
    // ppc->tvp_vac_num->print_all();


    // if tv-dev-icu-frac is not set, use default
    if (ppc->tvp_dev_icu_frac->empty_data()){
        ppc->tvp_dev_icu_frac->append_data(1.0);
    }
    if ( ppc->tvp_dev_icu_frac->fill_missing_days() ){
        ppc->tvp_dev_icu_frac->advance_indices(); // start day_index and data_index at 0
    } else {
        fprintf(stderr,"\n\n\tError in parsing time-varying dev-icu-frac.\n\tPlease check all -tv-dev-icu-frac-* command-line options or unset all these arguments to use the default values.\n\n");
        ppc->tvp_dev_icu_frac->print_all();
        exit(-1);
    }

    // if tv-dev-len-hospstay is not set, use default
    if (ppc->tvp_dev_len_hospstay->empty_data()){
        ppc->tvp_dev_len_hospstay->append_data(1.0);
    }
    if ( ppc->tvp_dev_len_hospstay->fill_missing_days() ){
        ppc->tvp_dev_len_hospstay->advance_indices(); // start day_index and data_index at 0
    } else {
        fprintf(stderr,"\n\n\tError in parsing time-varying dev-len-hospstay.\n\tPlease check all -tv-dev-len-hospstay-* command-line options or unset all these arguments to use the default values.\n\n");
        ppc->tvp_dev_len_hospstay->print_all();
        exit(-1);
    }
    // ppc->tvp_dev_len_hospstay->print_all();
    // exit(0);
    
    InitializeContactMatrices();
    
    
    /*for(int i=0; i < ppc->v_betas.size(); i++)
    {
        printf("\n\t time-start: %1.3f   \t beta-val: %1.3f", ppc->v_betatimes[i], ppc->v_betas[i] );
    }
    printf("\n\t time final is %1.1f\n", G_CLO_TF);
    
   return 1;*/
    
    
    
    //
    // ###  2.  INITIALIZE PARAMETERS - these are the default/starting values
    //
    ppc->v[ i_len_incub_period ]                            = 6.0;  // this is an average of the Lauer et al estimate (5.5d) and Backer et al (6.5d)
    ppc->v[ i_len_symptomatic_infectious_period_phase_1 ]   = G_CLO_TIME_SYMPTOHOSP;
    ppc->v[ i_len_symptomatic_infectious_period_phase_2 ]   = 7.0;
    // ppc->v[ i_len_medicalfloor_hospital_stay ]              = 10.7 * G_CLO_DEV_LEN_HOSPSTAY; //   take from Lewnard et al, survivors only
    if (ppc->tvp_dev_len_hospstay->get_day_index() == 0 && ppc->tvp_dev_len_hospstay->get_current_begin_day() > 0){
        // do nothing since the default dev-len-hospstay is 1.0
    } else
    {
        ppc->v[ i_len_medicalfloor_hospital_stay ] = 10.7 * ppc->tvp_dev_len_hospstay->get_current_data();
    }
    ppc->v[ i_mean_time_vent ]                              = G_CLO_MEANTIME_ON_VENT_SURV;   //   mean time on ventilator for survivors
    
    
    // params below are for relative infectiousness of certain individuals; I1 and I2 individuals have infectiousness = 1.0
    ppc->v[ i_phi_asymp ]           = 0.5;  // set to same value as Imperial college models, and JosephWu NatMed paper assumption (sensitivity analysis performed)
    ppc->v[ i_phi_symp_phase2 ]     = 0.75; // the second half of the symptomatic period is somewhat less infectious (TODO get evidence here)
    ppc->v[ i_phi_incub ]           = 0.5;  // currently a complete unknown
    ppc->v[ i_phi_hosp ]            = 1.0;
    ppc->v[ i_phi_hosp_recovering ] = 0.2;  // likely to be at late stage of infection with much less culturable virus
    ppc->v[ i_phi_icu ]             = 1.0;
    ppc->v[ i_phi_vent ]            = 0.1;  // to discuss
    
    
    ppc->v[ i_selfisolation_factor ] = G_CLO_SELFISOLATION_FACTOR;  // set to 1.0 by default; this means that infected and symptomatic individuals in the community 
                                                                    // do not self-isolate any more than anyone else; set this to something smaller than one
                                                                    // to make these individuals isolate more than non-symp people
    
    
    // params below are for relative contact levels of hospitalized, ICUed, and vented individuals
    
    int num_days_to_average=10;
    double betasum=0.0;
    for( int d=1; d<=num_days_to_average; d++) betasum += ppc->v_betas[d]; // ignore day 0 because this is the pre-March 1 period
    ppc->v[ i_prelockdown_beta ] = betasum / ((double) num_days_to_average );
    
    ppc->v[ i_beta_hosp ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v[ i_prelockdown_beta ];  // these relbeta's do not change with social distancing measures
    ppc->v[ i_beta_icu ]  = G_CLO_RELATIVE_BETA_HOSP * ppc->v[ i_prelockdown_beta ];  // because the contact rate of a hospitalized patient is unaffected by the 
    ppc->v[ i_beta_vent ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v[ i_prelockdown_beta ];  // social distancing policies set for healthy individuals
    
    ppc->v[ i_min_relbeta_00 ] = G_CLO_MIN_MIXINGLEVEL_00; // this is the minimum allowable relative beta for the 0-9 age group
    ppc->v[ i_min_relbeta_10 ] = G_CLO_MIN_MIXINGLEVEL_10; // this is the minimum allowable relative beta for the 10-19 age group
    ppc->v[ i_min_relbeta_20 ] = G_CLO_MIN_MIXINGLEVEL_20; // this is the minimum allowable relative beta for the 20-29 age group
    ppc->v[ i_min_relbeta_30 ] = G_CLO_MIN_MIXINGLEVEL_30; // this is the minimum allowable relative beta for the 30-39 age group
    ppc->v[ i_min_relbeta_40 ] = G_CLO_MIN_MIXINGLEVEL_40; // this is the minimum allowable relative beta for the 40-49 age group
    ppc->v[ i_min_relbeta_50 ] = G_CLO_MIN_MIXINGLEVEL_50; // this is the minimum allowable relative beta for the 50-59 age group
    ppc->v[ i_min_relbeta_60 ] = G_CLO_MIN_MIXINGLEVEL_60; // this is the minimum allowable relative beta for the 60-69 age group
    ppc->v[ i_min_relbeta_70 ] = G_CLO_MIN_MIXINGLEVEL_70; // this is the minimum allowable relative beta for the 70-79 age group
    ppc->v[ i_min_relbeta_80 ] = G_CLO_MIN_MIXINGLEVEL_80; // this is the minimum allowable relative beta for the 80-89 age group
    
    
    
    // set the relative susceptibilities of the different age groups
    //for(int ac=0;ac<NUMAC;ac++) 
    //{
    ppc->v_rel_susc[0] = G_CLO_RELSUSC_0;
    ppc->v_rel_susc[1] = G_CLO_RELSUSC_10;
    ppc->v_rel_susc[2] = G_CLO_RELSUSC_20;
    ppc->v_rel_susc[3] = G_CLO_RELSUSC_30;
    ppc->v_rel_susc[4] = G_CLO_RELSUSC_40;
    ppc->v_rel_susc[5] = G_CLO_RELSUSC_50;
    ppc->v_rel_susc[6] = G_CLO_RELSUSC_60;
    ppc->v_rel_susc[7] = G_CLO_RELSUSC_70;
    ppc->v_rel_susc[8] = G_CLO_RELSUSC_80;
    //}

    // set the fraction of individuals who progress to asymptomatic infection
    double a=G_CLO_SYMP_FRAC;  //NOTE this number has to be somewhere between 0.1 and 0.3, closer to 0.3 probably; default is 0.25;
    //
//     ppc->v_prob_E_A[0] = 1.0 - a * 0.05; // 98.5% to 99.5% asymp       // these relative params are all taken from Joseph Wu et al, Nat Med, 2020
//     ppc->v_prob_E_A[1] = 1.0 - a * 0.08; // 97.6% to 99.2% asymp
//     ppc->v_prob_E_A[2] = 1.0 - a * 0.41; // 87.7% to 95.9% asymp
//     ppc->v_prob_E_A[3] = 1.0 - a * 1.00; // 70% to 90% asymp           // 30-39, the reference group
//     ppc->v_prob_E_A[4] = 1.0 - a * 1.35; // 59.5% to 86.5% asymp
//     ppc->v_prob_E_A[5] = 1.0 - a * 1.99; // 40.3% to 80.1% asymp 
//     ppc->v_prob_E_A[6] = 1.0 - a * 2.85; // 14.5% to 71.5% asymp 
//     ppc->v_prob_E_A[7] = 1.0 - a * 3.05; //  8.5% to 69.5% asymp 
//     ppc->v_prob_E_A[8] = 1.0 - a * 2.49; // 25.3% to 75.1% asymp 

//     ppc->v_prob_E_A[0] = 0.83;      // these values are taken from a range of studies, and approximated 
//     ppc->v_prob_E_A[1] = 0.83;      // their current state (2020-06-09) is for a RI manual fit exercise
//     ppc->v_prob_E_A[2] = 0.5; 
//     ppc->v_prob_E_A[3] = 0.5; 
//     ppc->v_prob_E_A[4] = 0.45; 
//     ppc->v_prob_E_A[5] = 0.35; 
//     ppc->v_prob_E_A[6] = 0.35; 
//     ppc->v_prob_E_A[7] = 0.25; 
//     ppc->v_prob_E_A[8] = 0.2; 
    
    
    if( G_B_SYMP_FRAC_SIMPLEAVERAGE )
    {
        ppc->v_prob_E_A[0] = 0.381; // 0.83;      // these values have been updated (2020/07/08) using Fuhan's most recent lit review
        ppc->v_prob_E_A[1] = 0.381; // 0.83;      // excluding studies that did not distinguish between pre-symp and asymp, and doing a simple weighted average by sample size
        ppc->v_prob_E_A[2] = 0.460; // 0.5;       // the age-specific (in 20-yr ge bands) asymp fractions are
        ppc->v_prob_E_A[3] = 0.460; // 0.5;       // asymp fractions are 0.381 (0-19yo), 0.460 (20-39yo), 0.431 (40-59yo), 0.409 (60-79yo), 0.130 (80+) 
        ppc->v_prob_E_A[4] = 0.431; // 0.45; 
        ppc->v_prob_E_A[5] = 0.431; // 0.35; 
        ppc->v_prob_E_A[6] = 0.409; // 0.35; 
        ppc->v_prob_E_A[7] = 0.409; // 0.25; 
        ppc->v_prob_E_A[8] = 0.130; // 0.20;         
    }
    else if( G_B_SYMP_FRAC_DAVIES_20200616 )
    {
        ppc->v_prob_E_A[0] = 0.71;      // These asymptomatic fractions are taken from the Davies et al Nature Medicine paper; https://www.nature.com/articles/s41591-020-0962-9
        ppc->v_prob_E_A[1] = 0.79;      // Extended Data Fig 4, mean values for clinical fractions
        ppc->v_prob_E_A[2] = 0.73; 
        ppc->v_prob_E_A[3] = 0.67; 
        ppc->v_prob_E_A[4] = 0.60; 
        ppc->v_prob_E_A[5] = 0.51; 
        ppc->v_prob_E_A[6] = 0.37; 
        ppc->v_prob_E_A[7] = 0.31;      // Davies et al had one estimate for the >70 age group 
        ppc->v_prob_E_A[8] = 0.31;         
        
    }
    else // this means we are using G_CLO_SYMP_FRAC_EQUAL for all ages
    {
        for(int ac=0;ac<NUMAC;ac++) ppc->v_prob_E_A[ac] = 1.0 - G_CLO_SYMP_FRAC_EQUAL;
    }
    
    
    
    
    
    

    // set the fraction of individuals who are hospitalized immediately after I_2
    //
//     double b1 = G_CLO_HOSPFRAC_YOUNG_DEV;   // default is 1.8 :: REASON is that we want these rates to match the hosp-age-dist in the Lewnard paper
//     double b2 = G_CLO_HOSPFRAC_MID_DEV;     // default is 
//     double b3 = G_CLO_HOSPFRAC_OLD_DEV;     // default is 
/*     ppc->v_prob_I2_H[0] = G_CLO_HOSPFRAC_10;    // NOTE because there is not likely to be data on hosp rates for 0-9 y.o., we adopt the hosp rate for 10-19 year olds here
    ppc->v_prob_I2_H[1] = G_CLO_HOSPFRAC_10;
    ppc->v_prob_I2_H[2] = G_CLO_HOSPFRAC_20;
    ppc->v_prob_I2_H[3] = G_CLO_HOSPFRAC_30;    // this is the hospitalization probability for 30-39 year-olds
    ppc->v_prob_I2_H[4] = G_CLO_HOSPFRAC_40;
    ppc->v_prob_I2_H[5] = G_CLO_HOSPFRAC_50;
    ppc->v_prob_I2_H[6] = G_CLO_HOSPFRAC_60;
    ppc->v_prob_I2_H[7] = G_CLO_HOSPFRAC_70;
    ppc->v_prob_I2_H[8] = G_CLO_HOSPFRAC_80; */

    
    
    
    //
    // IMPORTANT-NOTE:  now that you have assigned the hospitalization probabilities, adjust them upwards for the 
    //                  early March period when hosp rates were higher than normal
    ppc->earlymarch_highhosp_factor = G_CLO_EARLYMARCH_HOSPRATE;
    ppc->earlymarch_highhosp_endday = G_CLO_EARLYMARCH_ENDDAY;
    //ppc->apply_earlymarch_hosprates();



    // set the probability of death for I4 for all 9 age classes
    ppc->v_prob_I4_D[0] = 0.0; // COMPLETELY UNKNOWN
    ppc->v_prob_I4_D[1] = 0.0;
    ppc->v_prob_I4_D[2] = 0.0;
    ppc->v_prob_I4_D[3] = 0.0;
    ppc->v_prob_I4_D[4] = 0.0;
    ppc->v_prob_I4_D[5] = 0.0;
    ppc->v_prob_I4_D[6] = G_CLO_DEATHPROB_HOME_60;
    ppc->v_prob_I4_D[7] = G_CLO_DEATHPROB_HOME_70;
    ppc->v_prob_I4_D[8] = G_CLO_DEATHPROB_HOME_80;

    // set the probability of progression from the "HA_1" -state to the CA state; this is the only way to go from hospitalizatio to ICU
    ppc->v_prob_HA_CA[0] = 0.304; 
    ppc->v_prob_HA_CA[1] = 0.293;     // NOTE the estimated Lewnard et al (medRxiv, Apr 16) probabilities are used here
    ppc->v_prob_HA_CA[2] = 0.2825;    //      averaged over male/female equally
    ppc->v_prob_HA_CA[3] = 0.301;
    ppc->v_prob_HA_CA[4] = 0.463;
    ppc->v_prob_HA_CA[5] = 0.4245;
    ppc->v_prob_HA_CA[6] = 0.460;
    ppc->v_prob_HA_CA[7] = 0.4835;
    ppc->v_prob_HA_CA[8] = 0.416;
    
    // the probabilities above range from: 
    // double c1 = G_CLO_ICUFRAC_DEV;
    for(int ac=0; ac<NUMAC; ac++)
    {
        // ppc->v_prob_HA_CA[ac] *= c1;
        //ppc->v_prob_HA_CA[ac] = 0.0;
        if (ppc->tvp_dev_icu_frac->get_day_index() == 0 && ppc->tvp_dev_icu_frac->get_current_begin_day() > 0){
            // do nothing since the default dev-icu-frac is 1.0
        } else
        {
            ppc->v_prob_HA_CA[ac] *= ppc->tvp_dev_icu_frac->get_current_data();
        }
        
    }
    
    // set the probability of death for HA4 for all 9 age classes - 
    ppc->v_prob_HA4_D[0] = 0.0;
    ppc->v_prob_HA4_D[1] = 0.0;
    ppc->v_prob_HA4_D[2] = 0.0;
    ppc->v_prob_HA4_D[3] = 0.0;
    ppc->v_prob_HA4_D[4] = 0.0;
    ppc->v_prob_HA4_D[5] = 0.0;   // NOTE there are no data right now for these numbers
    ppc->v_prob_HA4_D[6] = 0.0;
    ppc->v_prob_HA4_D[7] = G_CLO_PROB_NONICU_DEATH_80 / 2.0;
    ppc->v_prob_HA4_D[8] = G_CLO_PROB_NONICU_DEATH_80; 

    
    // set the probability of ventilation for CA-individuals for all 9 age classes - 
    ppc->v_prob_CA_V[0] = G_CLO_PROB_ICU_TO_VENT; // these can be set to about .75 for the higher age classes; from the Seattle ICU paper on 24 patients
    ppc->v_prob_CA_V[1] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[2] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[3] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[4] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[5] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[6] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[7] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[8] = G_CLO_PROB_ICU_TO_VENT;

    // from the Seattle ICU data on 24 patients, this probability is 60% -- obviously, it's a small sample size of older patients
    // these are being set to the ICU-to-Death probabilities since it's very difficult to get good data on death when on and not on a ventilator (for ICU patients)
    //
    double vd = G_CLO_VENTDEATH_MID_DEV; // default set to 0.7
    ppc->v_prob_V_D[0] = 0.03125;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[1] = 0.05119;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[2] = 0.15;          // this range should be between 14% (Lewnard) and 16.7%  (Graselli)
    ppc->v_prob_V_D[3] = 0.15;          // this range should be between 13% (Lewnard) and 17%  (Graselli), but Yang LRM observed 0%
    ppc->v_prob_V_D[4] = 0.400*vd;          // this range should be between 31% and 50%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[5] = 0.460*vd;          // this range should be between 26% and 70%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[6] = 0.585*vd;         // this range should be between 39% and 72%  (Lewnard, Graselli, Yang LRM, Bhatraju) 
    ppc->v_prob_V_D[7] = 0.70*G_CLO_VENTDEATH_70_DEV;          // this range should be between 60% and 88%  (Lewnard, Graselli, Yang LRM, Bhatraju)
    ppc->v_prob_V_D[8] = 0.90*G_CLO_VENTDEATH_80_DEV;          // this range should be between 60% and 100% (Lewnard, Graselli, Yang LRM, Bhatraju)
    

    // set the probability of death for CA-individuals for all 9 age classes 
    //
    // if you are in the ICU, and you don't progress to death, and your don't progress to ventilation, that means you are back on the medical-floor level of care
    //
    // very little data here; Seattle study on 24 ICU patients says this should be about 0.125
    //
    ppc->v_prob_CA_D[0] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];         
    ppc->v_prob_CA_D[1] = (1.0-ppc->v_prob_CA_V[1]) * ppc->v_prob_V_D[1];                                     
    ppc->v_prob_CA_D[2] = (1.0-ppc->v_prob_CA_V[2]) * ppc->v_prob_V_D[2];
    ppc->v_prob_CA_D[3] = (1.0-ppc->v_prob_CA_V[3]) * ppc->v_prob_V_D[3];
    ppc->v_prob_CA_D[4] = (1.0-ppc->v_prob_CA_V[4]) * ppc->v_prob_V_D[4];
    ppc->v_prob_CA_D[5] = (1.0-ppc->v_prob_CA_V[5]) * ppc->v_prob_V_D[5];
    ppc->v_prob_CA_D[6] = (1.0-ppc->v_prob_CA_V[6]) * ppc->v_prob_V_D[6];  // 0.146 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[7] = (1.0-ppc->v_prob_CA_V[7]) * ppc->v_prob_V_D[7];  // 0.175 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[8] = (1.0-ppc->v_prob_CA_V[8]) * ppc->v_prob_V_D[8];  // 0.225 --- so it's not close to the average of 0.125 obersved in Seattle 24-patent study; but these patients are >80
    
    
    // PROBABILITY OF DYING AFTER EXTUBATION, BUT WHILE STILL IN ICU
    ppc->v_prob_CR_D[0] = 0.0; 
    ppc->v_prob_CR_D[1] = 0.0;
    ppc->v_prob_CR_D[2] = 0.0;
    ppc->v_prob_CR_D[3] = 0.0;
    ppc->v_prob_CR_D[4] = 0.0;
    ppc->v_prob_CR_D[5] = 0.0;
    ppc->v_prob_CR_D[6] = 0.5*G_CLO_DEATHPROB_POSTVENT;     // this probability is 0.189 in the Gupta et al JAMA paper, mostly for older individuals
    ppc->v_prob_CR_D[7] = 1.0*G_CLO_DEATHPROB_POSTVENT;
    ppc->v_prob_CR_D[8] = 2.0*G_CLO_DEATHPROB_POSTVENT;
    
    
    
    
    
    //
    // ###  3.  RUN THE MODEL
    //
    if( !G_B_DEATHRATE_OUTPUT )
    {
        generate_trajectories( 0.01, 3.0, 0.0, G_CLO_TF, 0.0 );
        if( OutFile != NULL ) fclose(OutFile);
    }



    
    
    //
    // ###  4.  OUTPUT DIAGNOSTICS
    //
    if( G_B_DIAGNOSTIC_MODE )
    {
        int ac;
        
        printf("\n\t\t\t\t 0-9 \t\t 10-19  \t 20-29  \t 30-39  \t 40-49  \t 50-59  \t 60-69  \t 70-79  \t 80+ ");
        printf("\n\tTotal Symp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTJ + ac]+0.5) );
        }
        printf("\n\tTotal Deaths    ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTD + ac]+0.5) + (int)(yic[STARTDHOSP + ac]+0.5));
        }
        printf("\n\tCFR             ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t %1.3f%%   ", 100.0 * (yic[STARTD + ac]+yic[STARTDHOSP + ac]) / yic[STARTJ + ac]  );
        }

        printf("\n\n\tTotal Hosp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTK + ac]+0.5) );
        }
        printf("\n\tHosp FR         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTDHOSP + ac] / yic[STARTK + ac]  );
        }

        double total_num_hospitalized=0.0;
        for(ac=0; ac<NUMAC; ac++) total_num_hospitalized += yic[STARTK + ac];
        
        printf("\n\n\t%% Hosp by Age         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTK + ac] / total_num_hospitalized  );
        }
        printf("\n\tLewnard et al  \t\t  0.1%% \t\t  0.2%% \t\t  3.6%% \t\t  9.8%% \t\t  14.6%% \t  20.5%% \t  22.2%% \t  16.6%% \t  12.3%%");
        
    
        printf("\n\n");

        /*for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
        {
            for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
            {
                printf("\t %1.2f ", G_C_COMM[acs][aci]);
            }
            printf("\n");
        }
        
        printf("\n\n");*/
        
    }
    
    
    
    //
    // ###  5.  OUTPUT DEATH RATES
    //
    if( G_B_DEATHRATE_OUTPUT )
    {
        // TODO!!!
        // fprintf(stderr,"\n\n\tNeed to take care of time-varying hosp-frac (replacing v_prob_I2_H)!!!\n\n");
        // exit(-1);

        // note: v_prob_HA_CA uses the first dev-icu-frac

        double HFR[NUMAC];
        double CFR[NUMAC];
        double IFR[NUMAC];
        
        int ac;
        
        for(ac=0; ac<NUMAC; ac++)
        {
            double prob_death_given_icu = ppc->v_prob_CA_D[ac]  +  ppc->v_prob_CA_V[ac]*ppc->v_prob_V_D[ac]  +  (1.0-ppc->v_prob_CA_V[ac]-ppc->v_prob_CA_D[ac])*ppc->v_prob_HA4_D[ac];
            HFR[ac] = ppc->v_prob_HA_CA[ac]*prob_death_given_icu  +  (1.0-ppc->v_prob_HA_CA[ac])*ppc->v_prob_HA4_D[ac];
            
            CFR[ac] = ppc->v_prob_I2_H[ac]*HFR[ac] + (1.0-ppc->v_prob_I2_H[ac])*ppc->v_prob_I4_D[ac];
            
            IFR[ac] = (1.0-ppc->v_prob_E_A[ac])*CFR[ac];
            
        }
        
        
        if( OutFile == NULL )
        {
            
            double allages_ifr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_ifr += (yic[ac] / ppc->v[i_N]) * IFR[ac];
            double allages_cfr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_cfr += (yic[ac] / ppc->v[i_N]) * CFR[ac];
            double allages_hfr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_hfr += (yic[ac] / ppc->v[i_N]) * HFR[ac];

            for(ac=0; ac<NUMAC; ac++) printf("%1.6f \t ", IFR[ac] );
            for(ac=0; ac<NUMAC; ac++) printf("%1.5f \t ", CFR[ac] );
            for(ac=0; ac<NUMAC; ac++) printf("%1.4f \t ", HFR[ac] );
        
            printf("%1.6f \t %1.5f \t %1.4f", allages_ifr, allages_cfr, allages_hfr );
        }
        else
        {
            double allages_ifr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_ifr += (yic[ac] / ppc->v[i_N]) * IFR[ac];
            double allages_cfr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_cfr += (yic[ac] / ppc->v[i_N]) * CFR[ac];
            double allages_hfr=0.0;
            for(ac=0; ac<NUMAC; ac++) allages_hfr += (yic[ac] / ppc->v[i_N]) * HFR[ac];

            for(ac=0; ac<NUMAC; ac++) fprintf(OutFile, "%1.6f \t ", IFR[ac] );
            for(ac=0; ac<NUMAC; ac++) fprintf(OutFile, "%1.5f \t ", CFR[ac] );
            for(ac=0; ac<NUMAC; ac++) fprintf(OutFile, "%1.4f \t ", HFR[ac] );
        
            fprintf(OutFile, "%1.6f \t %1.5f \t %1.4f", allages_ifr, allages_cfr, allages_hfr );
            
            fclose(OutFile);
            
        }
        
        
    }
   

    delete[] yic;
    delete[] s_to_r;
    delete ppc;
    return 0;
}



void InitializeContactMatrices( void )
{
    // ### GENERAL POPULATION CONTACT MATRIX
    // if using symmetric contact matrix
    if ( !ppc->tvp_contact_rate->empty_data() ){
        // if the first tv-contact-rate-beginday CLO is higher than 0, use the default contact rate (1.0 for all groups) first 
        bool using_default_contact_rate_first = (ppc->tvp_contact_rate->get_day_index() == 0) && (ppc->tvp_contact_rate->get_current_begin_day() > 0);
        /* if (using_default_contact_rate_first) {printf("using_default_contact_rate_first\n"); }
        else {
            printf("#### NOT using_default_contact_rate_first ####\n");
            ppc->tvp_contact_rate->print_current_state();
        } */
        
        for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
        {
            for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
            {
                if (using_default_contact_rate_first){
                    G_C_COMM[acs][aci] = 1.0;
                } else {
                    G_C_COMM[acs][aci] = ppc->tvp_contact_rate->get_current_data(acs) * ppc->tvp_contact_rate->get_current_data(aci);
                }
                // G_C_COMM[acs][aci] = ppc->v_mixing_level[acs] * ppc->v_mixing_level[aci];
            }
        }
    } else if ( !ppc->tvp_contact_coeff->empty_data() )
    {
        // CoMix - Belgium wave 1 (late April 2020) matrix from `socialmixr` R package;
        G_C_COMM[0][0] = 0.1;

        G_C_COMM[0][1] = 0.161290322580645;     G_C_COMM[1][0] = 0.161290322580645; 
        G_C_COMM[0][2] = 0.272881355932203;     G_C_COMM[2][0] = 0.272881355932203;
        G_C_COMM[0][3] = 1.04761904761905;      G_C_COMM[3][0] = 1.04761904761905;
        G_C_COMM[0][4] = 0.472573839662447;     G_C_COMM[4][0] = 0.472573839662447;
        G_C_COMM[0][5] = 0.118043844856661;     G_C_COMM[5][0] = 0.118043844856661;
        G_C_COMM[0][6] = 0.169491525423729;     G_C_COMM[6][0] = 0.169491525423729;
        G_C_COMM[0][7] = 0.064814814814815;     G_C_COMM[7][0] = 0.064814814814815;
        G_C_COMM[0][8] = 0.001;                 G_C_COMM[8][0] = 0.001;

        G_C_COMM[1][1] = 1.95967741935484; 
        G_C_COMM[1][2] = 1.20967741935484;      G_C_COMM[2][1] = 0.389830508474576;
        G_C_COMM[1][3] = 0.564516129032258;     G_C_COMM[3][1] = 0.328407224958949;
        G_C_COMM[1][4] = 0.983870967741935;     G_C_COMM[4][1] = 1.21097046413502;
        G_C_COMM[1][5] = 0.491935483870968;     G_C_COMM[5][1] = 0.450252951096121;
        G_C_COMM[1][6] = 0.080645161290323;     G_C_COMM[6][1] = 0.149152542372881;
        G_C_COMM[1][7] = 0.258064516129032;     G_C_COMM[7][1] = 0.037037037037037;
        G_C_COMM[1][8] = 0.001;                 G_C_COMM[8][1] = 0.001;

        G_C_COMM[2][2] = 1.93728813559322;
        G_C_COMM[2][3] = 0.377966101694915;     G_C_COMM[3][2] = 1.1559934318555;
        G_C_COMM[2][4] = 0.961016949152542;     G_C_COMM[4][2] = 0.727144866385373;
        G_C_COMM[2][5] = 0.8;                   G_C_COMM[5][2] = 0.860033726812816;
        G_C_COMM[2][6] = 0.111864406779661;     G_C_COMM[6][2] = 0.688135593220339;
        G_C_COMM[2][7] = 0.133898305084746;     G_C_COMM[7][2] = 0.236111111111111;
        G_C_COMM[2][8] = 0.150847457627119;     G_C_COMM[8][2] = 0.001;

        G_C_COMM[3][3] = 1.26929392446634;
        G_C_COMM[3][4] = 0.48111658456486;      G_C_COMM[4][3] = 1.0365682137834;
        G_C_COMM[3][5] = 0.41871921182266;      G_C_COMM[5][3] = 0.458684654300169;
        G_C_COMM[3][6] = 0.183908045977011;     G_C_COMM[6][3] = 0.603389830508475;
        G_C_COMM[3][7] = 0.139573070607553;     G_C_COMM[7][3] = 0.324074074074074;
        G_C_COMM[3][8] = 0.004926108374384;     G_C_COMM[8][3] = 0.153846153846154;

        G_C_COMM[4][4] = 1.06188466947961;
        G_C_COMM[4][5] = 0.312236286919831;     G_C_COMM[5][4] = 0.728499156829679;
        G_C_COMM[4][6] = 0.135021097046413;     G_C_COMM[6][4] = 0.415254237288136;
        G_C_COMM[4][7] = 0.174402250351617;     G_C_COMM[7][4] = 0.449074074074074;
        G_C_COMM[4][8] = 0.035161744022504;     G_C_COMM[8][4] = 0.692307692307692;

        G_C_COMM[5][5] = 0.772344013490725;
        G_C_COMM[5][6] = 0.124789207419899;     G_C_COMM[6][5] = 0.844067796610169;
        G_C_COMM[5][7] = 0.12141652613828;      G_C_COMM[7][5] = 0.37962962962963;
        G_C_COMM[5][8] = 0.084317032040472;     G_C_COMM[8][5] = 0.153846153846154;

        G_C_COMM[6][6] = 0.691525423728813;
        G_C_COMM[6][7] = 0.347457627118644;     G_C_COMM[7][6] = 0.574074074074074;
        G_C_COMM[6][8] = 0.030508474576271;     G_C_COMM[8][6] = 0.230769230769231;

        G_C_COMM[7][7] = 1.40740740740741;
        G_C_COMM[7][8] = 0.083333333333333;     G_C_COMM[8][7] = 0.307692307692308;

        G_C_COMM[8][8] = 0.384615384615385;

        bool using_default_contact_coeff_first = (ppc->tvp_contact_coeff->get_day_index() == 0) && (ppc->tvp_contact_coeff->get_current_begin_day() > 0);

        if ( !using_default_contact_coeff_first ){
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
        }
    } else
    {
        fprintf(stderr,"\n\n\tNot enough information to construct contact matrices.\n\n");
        exit(-1);
    }
    
    

    
    // ### CONTACT MATRIX BTW INFECTED HOSPITALIZED INDIVIDUALS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs==0 )        // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.0; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; HEALTHY CHILDREN DO NOT VISIT THE HOSPITAL
            }
            else if( acs==1 )   // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.0; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; HEALTHY TEENS VISIT THE HOSPITAL RARELY
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.2; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; SOME MAY WORK IN THE HOSPITAL
            }
            else
            {
                G_C_HOSP[acs][aci] = 1.0;
            }
        }
    }

    
    // ### CONTACT MATRIX BTW INFECTED ICU PATIENTS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs<=1 )        // age-class of the susceptible contact
            {
                G_C_ICU[acs][aci] = 0.0; // HEALTHY CHILDREN & TEENS DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_ICU[acs][aci] = 0.0; // UNINFECTED ELDERLY DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else
            {
                G_C_ICU[acs][aci] = 1.0;
            }
        }
    }
    
    
    // ### CONTACT MATRIX BTW INFECTED-AND-VENTILATED ICU PATIENTS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs<=1 )        // age-class of the susceptible contact
            {
                G_C_VENT[acs][aci] = 0.0; // HEALTHY CHILDREN & TEENS DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_VENT[acs][aci] = 0.0; // UNINFECTED ELDERLY DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else
            {
                G_C_VENT[acs][aci] = 1.0;
            }
        }
    }
    
    
    
}


