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
#include "time_varying_prms_nonage.h"

#include <ctime>
#include <sys/time.h>

using namespace std;


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double b,d,s,c0,c2,c3;	// this are defined, i believe, in rkf.cpp
double* yic;

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

double G_CLO_MIXINGLEVEL[NUMAC];
double G_CLO_MIXINGLEVEL_POSTLD[NUMAC];     // the age-specific mixing level post-lockdown
int G_CLO_FIRSTLOCKDOWN_ENDDAY = 121;       // this is April 30 2020
bool G_CLO_POSTLD_MIXING_SET = false;

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

double G_C_COMM[NUMAC][NUMAC]; // this is the contact rate in the community between susc (first index) and infected (second index)                              - this matrix is symmetric
double G_C_HOSP[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected hospitalized patient (second index)     - this matrix is not symmetric
double G_C_ICU[NUMAC][NUMAC];  // this is the contact rate btw the community (first index, susceptible) and an infected ICU patient (second index)              - this matrix is not symmetric
double G_C_VENT[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected Ventilated patient (second index)       - this matrix is not symmetric
 

// END ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void InitializeContactMatrices( void );


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///////////////////// main flow ////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    time_varying_prms_nonage<double> *tvp_nonage = new time_varying_prms_nonage<double>();
    
    //
    // ###  1.  ALLOCATE SPACE FOR A PARAMETERS CLASS AND FOR THE MAIN STATE VARIABLES
    //
    
    ppc = new prms; assert( ppc );
    yic = new double[STARTK+NUMAC];
    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out

    // assign some default values to the v_beta and v_betatimes arrays
    ppc->v_betas.push_back( 1.0 );
    ppc->v_betas.push_back( 0.9 );
    ppc->v_betas.push_back( 0.8 );
    ppc->v_betatimes.push_back( 0.0 );
    ppc->v_betatimes.push_back( 60.0 );
    ppc->v_betatimes.push_back( 95.0 );
    


    // ppc->tvp_contact_rate->append_end_day(10);
    // ppc->tvp_contact_rate->append_end_day(20);
    // for (size_t i = 0; i < NUMAC; i++){ ppc->tvp_contact_rate->append_data(i, i*1.5); }
    



    // INITIALIZE GLOBAL VARIABLES
    
    // ppc->tvp_contact_rate->reset();

    ParseArgs( argc, argv );

    tvp_nonage->append_data(1213);
    tvp_nonage->fill_missing_days();
    tvp_nonage->print_all();
    tvp_nonage->advance_indices();
    printf("%d\n", tvp_nonage->get_data_index() );


    // if no contact rate is set via command-line options, use default
    if (ppc->tvp_contact_rate->empty_data() && ppc->tvp_contact_coeff->empty_data() ) {
        for (size_t i = 0; i < NUMAC; i++) { ppc->tvp_contact_rate->append_data(i, 1.0); }
    } else if (!ppc->tvp_contact_rate->empty_data() && !ppc->tvp_contact_coeff->empty_data() ) {
        fprintf(stderr,"\n\n\t-tv-contact-rate-* and -tv-contact-coeff-* command-line options are mutually exclusive.\n\n");
        ppc->tvp_contact_rate->print_all();
        ppc->tvp_contact_coeff->print_all();
    }
    
    // fill missing begin and/or end days in tvp_contact_rate
    // advance if good_days() and good_integrity()
    // else, print error message
    if ( ppc->tvp_contact_rate->fill_missing_days() ){
        ppc->tvp_contact_rate->advance_indices(); // start day_index and data_index at 0
    } else {
        fprintf(stderr,"\n\n\tError in parsing contact rates.\n\tPlease check all -tv-contact-rate-* command-line options or unset all these arguments to use the default values.\n\n");
        ppc->tvp_contact_rate->print_all();
    }
    

    printf("check integrity tvp_contact_rate %d\n", ppc->tvp_contact_rate->good_integrity() );

    ppc->tvp_contact_rate->print_all();
    // ppc->tvp_contact_rate->print_current_state();
    // ppc->tvp_contact_rate->advance_indices();
    // ppc->tvp_contact_rate->print_current_state();
    // ppc->tvp_contact_rate->reset();
    // ppc->tvp_contact_rate->print_all();


    // if no contact rate is set via command-line options, use default
    if (ppc->tvp_hosp_frac->empty_data()) {
        for (size_t i = 0; i < NUMAC; i++) { ppc->tvp_hosp_frac->append_data(i, 1.0); }
    }
    // fill missing begin and/or end days in tvp_hosp_frac
    // advance if good_days() and good_integrity()
    // else, print error message
    if ( ppc->tvp_hosp_frac->fill_missing_days() ){
        ppc->tvp_hosp_frac->advance_indices(); // start day_index and data_index at 0
    } else {
        fprintf(stderr,"\n\n\tError in parsing hospitalization fractions.\n\tPlease check all -tv-hosp-frac-* command-line options or unset all these arguments to use the default values.\n\n");
        ppc->tvp_hosp_frac->print_all();
    }
    
    printf("check integrity tvp_hosp_frac %d\n", ppc->tvp_hosp_frac->good_integrity() );

    ppc->tvp_hosp_frac->print_all();


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
            fprintf(stderr,"\n\n\tError in parsing contact rates.\n\tPlease check all -tv-vaccinees-* command-line options or unset all these arguments to use the default values.\n\n");
            ppc->tvp_vac_num->print_all();
        }
    }
    printf("\ntvp_vac_num\n");
    ppc->tvp_vac_num->print_all();
        
    InitializeContactMatrices();
    
    

    delete[] yic;
    delete ppc;
    return 0;
}



void InitializeContactMatrices( void )
{
    // ### GENERAL POPULATION CONTACT MATRIX
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