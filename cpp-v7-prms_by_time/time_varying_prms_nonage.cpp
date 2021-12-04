#include <assert.h>
#include <stdio.h>
// #include "time_varying_prms.h"
#include "time_varying_prms_nonage.h"

extern double G_CLO_TF;

template <typename T>
time_varying_prms_nonage<T>::time_varying_prms_nonage() : time_varying_prms<T>() {
    v_data.clear();
    assert( v_data.size()==0 );
}

template <typename T>
time_varying_prms_nonage<T>::~time_varying_prms_nonage(){}

template <typename T>
void time_varying_prms_nonage<T>::reset(){
    this->v_begin_days.clear();
    this->v_end_days.clear();
    v_data.clear() ;
    this->day_index = -1;
    this->data_index = -1;

    assert( this->v_begin_days.size()==0 );
    assert( this->v_end_days.size()==0 );
    assert( v_data.size()==0 );
}


//// getters
template <typename T>
T time_varying_prms_nonage<T>::get_current_data() const{ 
    assert(this->data_index > -1 && this->data_index < (int)v_data.size()); 
    return v_data[this->data_index];
}
template <typename T>
T time_varying_prms_nonage<T>::get_current_data(const size_t& age_group) const{ 
    fprintf(stderr,"\n\n\tThe data vector is not age-stratified.\n\n"); return -1;
}

//// setters
template <typename T>
void time_varying_prms_nonage<T>::append_data(const T& data){ v_data.emplace_back(data); }
template <typename T>
void time_varying_prms_nonage<T>::append_data(const size_t& age_group, const T& data){ fprintf(stderr,"\n\n\tThe data vector is not age-stratified.\n\n"); }


//// among v_begin_days and v_end_days, fill the shorter one with reasonable data
//// add the last pair of begin-end day if the last input endday is < G_CLO_TF
//// if both v_begin_days and v_end_days are empty but v_data is not, add 0 to v_begin_day and G_CLO_TF+1 to v_end_days
//// return true if both vectors are not empty and their values are good 
template <typename T>
bool time_varying_prms_nonage<T>::fill_missing_days(){
    if (empty_data()){
        return this->good_days() && good_integrity();
    }

    if (this->v_begin_days.empty() && this->v_end_days.empty() && !empty_data() ){
        this->v_begin_days.emplace_back( 0 );
        this->v_end_days.emplace_back( ((unsigned)G_CLO_TF) + 1 );
    } else if (this->v_begin_days.size() > this->v_end_days.size()){
        for (size_t i = this->v_end_days.size(); i < this->v_begin_days.size()-1; i++){
            this->v_end_days.emplace_back( this->v_begin_days[i+1]-1 );
        }
        this->v_end_days.emplace_back( ((unsigned)G_CLO_TF) + 1 );
    } else if (this->v_end_days.size() > this->v_begin_days.size()){
        for (size_t i = this->v_begin_days.size(); i < this->v_end_days.size(); i++){
            if (i==0){ this->v_begin_days.emplace_back( 0 ); } // case v_begin_days is empty, add 0 first
            else { this->v_begin_days.emplace_back( this->v_end_days[i-1]+1 ); }
        }
    }

    if (this->v_end_days.back() < G_CLO_TF){
        this->v_begin_days.emplace_back( this->v_end_days.back() + 1 );
        this->v_end_days.emplace_back( ((unsigned)G_CLO_TF) + 1 );
    }
    // if ( !good_days() ){printf("check days\n"); }
    // if ( !good_integrity() ){ printf("check integrity\n"); }
    return this->good_days() && good_integrity();
}


//// advance indices if possible
//// return false if either of the three core vectors are empty or both of the new indicies point to the last element in those vectors
////    i.e. calling advance_indices() again won't have any effect on either of the indices.
//// return true if either of the new indicies point to elements before the last in the three core vectors
////    i.e. a future advance_indices() call is meaningful.
template <typename T>
bool time_varying_prms_nonage<T>::advance_indices(){
    if (this->v_begin_days.empty() || this->v_end_days.empty() || empty_data() ){ 
        return false;
    } else {
        this->day_index += (short)(this->v_begin_days.size() > this->day_index+1);
        this->data_index += (short)(v_data.size() > this->data_index+1);
        assert( this->data_index <= this->day_index );
    }    
    return (this->v_begin_days.size() > this->day_index+1 || v_data.size() > this->data_index+1);
}


//// check empty vectors
template <typename T>
bool time_varying_prms_nonage<T>::empty_data(const size_t& age_group) const{ fprintf(stderr,"\n\n\tThe data vector is not age-stratified.\n\n");  return NULL ; }
template <typename T>
bool time_varying_prms_nonage<T>::empty_data() const{ return v_data.empty(); }


//// check integrity of internal data
//// v_begin_days and v_end_days must have the same length
//// nested vectors in v_data must have identical length
//// all indices must be smaller than the size of corresponding vector(s)
template <typename T>
bool time_varying_prms_nonage<T>::good_integrity() const{
    if ( empty_data() ){
        return (this->day_index < 0 && this->data_index < 0);
    }
    bool good = (this->v_begin_days.size() == this->v_end_days.size()) && (this->day_index < (int)this->v_begin_days.size());
    good = good && ( (this->v_begin_days.empty() || this->v_end_days.empty()) ? v_data.empty() : v_data.size() > 0);
    good = good && this->data_index < (int)v_data.size();
    
    // print_all();
    return good;
}


//// print all vector
template <typename T>
void time_varying_prms_nonage<T>::print_all(){
    printf("v_begin_days\t");
    for (size_t i = 0; i < this->v_begin_days.size(); i++){ printf("%d ", this->v_begin_days[i]); }
    printf("\nv_end_days\t");
    for (size_t i = 0; i < this->v_end_days.size(); i++){ printf("%d ", this->v_end_days[i]); }
    printf("\nv_data\t");
    for (size_t i = 0; i < v_data.size(); i++){ printf("%.5f ", v_data[i]); }
    printf("\n");
}

template <typename T>
void time_varying_prms_nonage<T>::print_current_state(){
    printf("v_begin_days\t%d\n", this->get_current_begin_day() );
    printf("v_end_days\t%d\n", this->get_current_end_day() );
    printf("v_data\t%.5f\n", get_current_data() );
}


//// force instantiating time_varying_prms_nonage<double>
template class time_varying_prms_nonage<double>;