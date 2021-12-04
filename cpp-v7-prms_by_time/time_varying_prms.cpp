#include <stdio.h>
#include <assert.h>
#include "essentials.h"
#include "time_varying_prms.h"

extern double G_CLO_TF;

// constructor
template<typename T>
time_varying_prms<T>::time_varying_prms(/* args */)
{
    v_begin_days.clear();
    v_end_days.clear();
    // v_data.insert( v_data.begin(), NUMAC, std::vector<double>{} );
    day_index = -1;
    data_index = -1;


    assert( v_begin_days.size()==0 );
    assert( v_end_days.size()==0 );
    // assert( v_data.size()==NUMAC );
    // for (size_t i = 0; i < NUMAC; i++){  assert(v_data[i].size()==0);   }
}

template<typename T>
void time_varying_prms<T>::reset(){
    v_begin_days.clear();
    v_end_days.clear();
    // for (size_t i = 0; i < NUMAC; i++){ v_data[i].clear() ;   }
    day_index = -1;
    data_index = -1;

    assert( v_begin_days.size()==0 );
    assert( v_end_days.size()==0 );
    // assert( v_data.size()==NUMAC );
    // for (size_t i = 0; i < NUMAC; i++){  assert(v_data[i].size()==0);   }
}

// destructor
template<typename T>
time_varying_prms<T>::~time_varying_prms()
{
    // std::vector<double>().swap( v_begin_days );
    // std::vector<double>().swap( v_end_days );
    // std::vector<std::vector<double>>().swap( v_data );
}

//// getters
template<typename T>
int time_varying_prms<T>::get_day_index() const{ return day_index; }
template<typename T>
int time_varying_prms<T>::get_data_index() const{ return data_index; }
template<typename T>
unsigned time_varying_prms<T>::get_current_begin_day() const{ return (day_index < 0) ? -1 : v_begin_days[day_index]; }
template<typename T>
unsigned time_varying_prms<T>::get_current_end_day() const{ return (day_index < 0) ? -1 : v_end_days[day_index]; }

//// setters
template<typename T>
void time_varying_prms<T>::append_begin_day(const unsigned& day){ v_begin_days.emplace_back(day); }
template<typename T>
void time_varying_prms<T>::append_end_day(const unsigned& day){ v_end_days.emplace_back(day); }


//// check empty vectors
template<typename T>
bool time_varying_prms<T>::empty_begin_days() const{ return v_begin_days.empty(); }
template<typename T>
bool time_varying_prms<T>::empty_end_days() const{ return v_end_days.empty(); }

//// check each pair of begin-end day
//// return true if v_begin_days and v_end_days are not empty and each begin_day must be smaller than the corresponding end_day
//// return false otherwise
template<typename T>
bool time_varying_prms<T>::good_days() const{
    assert( v_begin_days.size() == v_end_days.size() );
    bool good = v_begin_days.empty() ? false : v_begin_days[0] < v_end_days[0] ;
    for (size_t i = 1; i < v_begin_days.size(); i++){
        good = good && v_begin_days[i-1] < v_begin_days[i] && v_begin_days[i] < v_end_days[i];
        if (!good){ break; }
    }
    return good;
}


//// force instantiating time_varying_prms<double>
template class time_varying_prms<double>;
