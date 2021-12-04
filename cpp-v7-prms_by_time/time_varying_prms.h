#ifndef TIME_VARYING_PRMS
#define TIME_VARYING_PRMS

#include <vector>

template<typename T>
class time_varying_prms
{
protected:
    std::vector<unsigned> v_begin_days;
    std::vector<unsigned> v_end_days;
    // std::vector<std::vector<double>> v_data;
    int day_index;
    int data_index;
public:
    time_varying_prms(/* args */);
    virtual ~time_varying_prms() = 0;

    int get_day_index() const;
    int get_data_index() const;
    unsigned get_current_begin_day() const;
    unsigned get_current_end_day() const;
    virtual T get_current_data() const = 0;
    virtual T get_current_data(const size_t& age_group) const = 0;

    void append_begin_day(const unsigned& day);
    void append_end_day(const unsigned& day);
    virtual void append_data(const T& data) = 0;
    virtual void append_data(const size_t& age_group, const T& data) = 0;
    
    virtual bool fill_missing_days() = 0;    
    virtual bool advance_indices() = 0;
    bool empty_begin_days() const;
    bool empty_end_days() const;
    virtual bool empty_data() const = 0;
    virtual bool empty_data(const size_t& age_group) const = 0;
    bool good_days() const;
    virtual bool good_integrity() const = 0;

    virtual void print_all() = 0;
    virtual void print_current_state() = 0;

    void reset();
};

#endif // TIME_VARYING_PRMS