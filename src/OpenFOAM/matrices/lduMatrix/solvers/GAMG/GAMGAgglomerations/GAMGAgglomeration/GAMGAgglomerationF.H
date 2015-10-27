#pragma once

#include "gpuAtomic.H"

namespace Foam
{

template<class Type>
struct GAMGAgglomerationRestrictFunctor
{
    Type* cf;

    GAMGAgglomerationRestrictFunctor
    (
        Type* _cf
    ):
        cf(_cf)
    {}

    template<class Tuple>
    __device__
    void operator()(const Tuple& t)
    {
        const label& id = thrust::get<0>(t);
        const Type& val = thrust::get<1>(t);

        atomicAdd(cf+id,val);
    }
};


template<class Type>
struct GAMGAgglomerationRestrictFaceFunctor
{
    Type* cf;

    GAMGAgglomerationRestrictFaceFunctor
    (
        Type* _cf
    ):
        cf(_cf)
    {}

    template<class Tuple>
    __device__
    void operator()(const Tuple& t)
    {
        const label& id = thrust::get<0>(t);

        if(id >= 0)
        {
            const Type& val = thrust::get<1>(t);
            atomicAdd(cf+id,val);
        }
    }
};

}
