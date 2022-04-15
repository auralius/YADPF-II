#include "COptVar.h"


COptVar::COptVar()
{
    _state_fn = NULL;
    _length   = 0;
    _max_val  = 0.0;
    _min_val  = 0.0;   
}

COptVar::~COptVar()
{
}

void COptVar::define(const optvarfntype f, int length)
{
    _state_fn = f;
    _length   = length;
    _max_val  = eval_at(_length - 1);
    _min_val  = eval_at(0);
}

inline double COptVar::eval_at(int index)
{
    if (_state_fn)
        return _state_fn(index);
    else
        return 0.0;
};
