/**
 * \file COptVar.h
 * \brief The class for the optimization variable: the state and input variables.
 */

#pragma once


#include <cstddef>


typedef double (*optvarfntype)(int);  /// Discretized state/input function call


class COptVar
{
public:
    COptVar();
    ~COptVar();

    void define(const optvarfntype f, int length);

    inline unsigned length() { return _length; };
    inline double max() { return _max_val; };
    inline double min() { return _min_val; };

    double eval_at(int index);
    
private:
    int _length;
    double _max_val;
    double _min_val;
    
    optvarfntype _state_fn;     /// Pointer to function call for discretizing the state/input variables
};

