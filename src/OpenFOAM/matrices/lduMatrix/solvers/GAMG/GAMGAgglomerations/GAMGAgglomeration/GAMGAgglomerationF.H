#pragma once

namespace Foam
{
namespace GAMG
{


template<class Type>
struct restrict
{
    const Type* ff;
    const label* sort;
    const Type zero;

    restrict
    (
        const Type* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort),
        zero(pTraits<Type>::zero)
    {}

    __host__ __device__
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
struct atomicRestrict
{
    Type* cf;
    const Type* ff;
    const label* fineToCoarse;

    atomicRestrict
    (
        Type* _cf,
        const Type* _ff,
        const label* _fineToCoarse
    ):
        cf(_cf),
        ff(_ff),
        fineToCoarse(_fineToCoarse)
    {}

    __device__
    void operator()(const label i)
    {
        #if __CUDA_ARCH__ >= 600
        atomicAdd(cf+fineToCoarse[i], ff[i]);
        #endif
    }
};


template<class Type>
struct atomicFaceRestrict
{
    Type* cf;
    const Type* ff;
    const label* fineToCoarse;

    atomicFaceRestrict
    (
        Type* _cf,
        const Type* _ff,
        const label* _fineToCoarse
    ):
        cf(_cf),
        ff(_ff),
        fineToCoarse(_fineToCoarse)
    {}

    __device__
    void operator()(const label ffacei)
    {
        #if __CUDA_ARCH__ >= 600
        label cFace = fineToCoarsePtr[ffacei];
        if (cFace >= 0)
        {
            atomicAdd(cfPtr+cFace, ffPtr[ffacei]);
        }
        #endif
    }
};


template<class Type>
struct prolong
{
    Type* ff;
    const Type* cf;
    const label* sort;
    const label* target;
    const label* targetStart;

    prolong
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

    __host__ __device__
    void operator()(const label& id)
    {
        Type val = cf[target[id]];
     
        for(label i = targetStart[id]; i < targetStart[id+1]; i++)
        {
            ff[sort[i]] = val;
        }
    }
};


}
}
