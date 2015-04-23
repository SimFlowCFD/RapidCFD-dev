#pragma once

namespace Foam
{

template<class Type>
struct GAMGAgglomerationRestrictFunctor
{
    const Type* ff;
    const label* sort;
    const Type zero;

    GAMGAgglomerationRestrictFunctor
    (
        const Type* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort),
        zero(pTraits<Type>::zero)
    {}

    __HOST____DEVICE__
    Type operator()(const label& start, const label& end)
    {
        Type out = zero;

        for(label i = start; i<end; i++)
        {
            out += ff[sort[i]];
        }

        return out;
    }
};

template<class Type>
struct GAMGAgglomerationProlongFunctor
{
    Type* ff;
    const Type* cf;
    const label* sort;
    const label* target;
    const label* targetStart;

    GAMGAgglomerationProlongFunctor
    (
        Type* _ff,
        const Type* _cf,
        const label* _sort,
        const label* _target,
        const label* _targetStart
    ):
        ff(_ff),
        cf(_cf),
        sort(_sort),
        target(_target),
        targetStart(_targetStart)
    {}

    __HOST____DEVICE__
    void operator()(const label& id)
    {
        Type val = cf[target[id]];
     
        for(label i = targetStart[id]; i < targetStart[id+1]; i++)
        {
            ff[sort[i]] = val;
        }
    }
};


template<class Type>
struct nonNegativeGAMGFunctor
{
    __HOST____DEVICE__
    bool operator()(const Type& x)
    {
        return x >= 0;
    }
};

}
