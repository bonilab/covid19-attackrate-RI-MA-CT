#ifndef TIME_VARYING_PRMS_NONAGE
#define TIME_VARYING_PRMS_NONAGE

#include <vector>
#include "time_varying_prms.h"


template <typename T>
class time_varying_prms_nonage : public time_varying_prms<T>
{
private:
    std::vector<T> v_data;
public:
    time_varying_prms_nonage(/* args */);
    ~time_varying_prms_nonage();

    T get_current_data() const ;
    T get_current_data(const size_t& age_group) const ;

    void append_data(const T& data) ;
    void append_data(const size_t& age_group, const T& data) ;

    bool fill_missing_days() ;
    bool advance_indices() ;

    bool empty_data() const ;
    bool empty_data(const size_t& age_group) const ;

    bool good_integrity() const ;

    void print_all() ;
    void print_current_state() ;

    void reset();
};

#endif // TIME_VARYING_PRMS_NONAGE