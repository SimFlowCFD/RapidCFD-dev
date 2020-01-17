#pragma once

namespace Foam
{
namespace GAMG
{


struct asymAgglomerate
{
    const scalar* uf;
    const scalar* lf;
    const bool* flip;
    const label* sort;

    asymAgglomerate
    (
        const scalar* _uf,
        const scalar* _lf,
        const bool* _flip,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        flip(_flip),
        sort(_sort)
    {}
    
    template<class Tuple,class Tuple2>
    __host__ __device__
    thrust::tuple<scalar,scalar> operator()(const Tuple& in,const Tuple2& t)
    {
        scalar uc = thrust::get<0>(in);
        scalar lc = thrust::get<1>(in); 
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            bool flipVal = flip[index];
            scalar ufVal = uf[index];
            scalar lfVal = lf[index];
            if( ! flipVal)
            {
                uc += ufVal;
                lc += lfVal;
            }
            else
            {
                uc += lfVal;
                lc += ufVal;
            }
        }

        return thrust::make_tuple(uc,lc);
    }

};


struct diagAsymAgglomerate
{
    const scalar* uf;
    const scalar* lf;
    const label* sort;

    diagAsymAgglomerate
    (
        const scalar* _uf,
        const scalar* _lf,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        sort(_sort)
    {}

    template<class Tuple>
    __host__ __device__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            out += uf[index] + lf[index];
        }

        return out;
    }
};


struct symAgglomerate
{
    const scalar* ff;
    const label* sort;

    symAgglomerate
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}
    
    template<class Tuple>
    __host__ __device__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += ff[sort[i]];
        }

        return out;
    }
};


struct diagSymAgglomerate
{
    const scalar* ff;
    const label* sort;

    diagSymAgglomerate
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}

    template<class Tuple>
    __host__ __device__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += 2*ff[sort[i]];
        }

        return out;
    }
};


struct symAtomicAgglomerate
{
    scalar* coarseDiag;
    scalar* coarseUpper;
    const scalar* fineUpper;
    const label* faceRestrictAddr;

    symAtomicAgglomerate
    (
        scalar* _coarseDiag,
        scalar* _coarseUpper,
        const scalar* _fineUpper,
        const label* _faceRestrictAddr
    ):
        coarseDiag(_coarseDiag),
        coarseUpper(_coarseUpper),
        fineUpper(_fineUpper),
        faceRestrictAddr(_faceRestrictAddr)
    {}

    template<class T>
    __device__
    inline void impl(T* coarseDiag, T* coarseUpper, const T* fineUpper, const label& fineFacei)
    {
        label cFace = faceRestrictAddr[fineFacei];
        T fine = fineUpper[fineFacei];
        if (cFace < 0)
        {
            label cell = -1 - cFace;
            atomicAdd(coarseDiag + cell, 2*fine);
        }
        else
        {
            atomicAdd(coarseUpper + cFace, fine);
        }
    }

    __device__
    inline void invoke(float* coarseDiag, float* coarseUpper, const float* fineUpper, const label& fineFacei)
    {
        impl(coarseDiag, coarseUpper, fineUpper, fineFacei);
    }

    __device__
    inline void invoke(double* coarseDiag, double* coarseUpper, const double* fineUpper, const label& fineFacei)
    {
        #if __CUDA_ARCH__ >= 600
        impl(coarseDiag, coarseUpper, fineUpper, fineFacei);
        #endif
    }

    __device__
    inline void operator()(const label& fineFacei)
    {
        invoke(coarseDiag, coarseUpper, fineUpper, fineFacei);
    }
};


struct asymAtomicAgglomerate
{
    scalar* coarseDiag;
    scalar* coarseUpper;
    scalar* coarseLower;
    const scalar* fineUpper;
    const scalar* fineLower;
    const label* faceRestrictAddr;
    const bool* faceFlip;

    asymAtomicAgglomerate
    (
        scalar* _coarseDiag,
        scalar* _coarseUpper,
        scalar* _coarseLower,
        const scalar* _fineUpper,
        const scalar* _fineLower,
        const label* _faceRestrictAddr,
        const bool* _faceFlip
    ):
        coarseDiag(_coarseDiag),
        coarseUpper(_coarseUpper),
        coarseLower(_coarseLower),
        fineUpper(_fineUpper),
        fineLower(_fineLower),
        faceRestrictAddr(_faceRestrictAddr),
        faceFlip(_faceFlip)
    {}

    template<class T>
    __device__
    inline void impl(T* cDiag, T* cUpper, T* cLower, const T* fUpper, const T* fLower, const label& facei)
    {
        label cFace = faceRestrictAddr[facei];
        T fu = fUpper[facei];
        T fl = fLower[facei];
        if (cFace < 0)
        {
            label cell = -1 - cFace;
            atomicAdd(coarseDiag + cell, fu + fl);
        }
        else
        {
            T* cu = cUpper + cFace;
            T* cl = cLower + cFace;
            bool flip = faceFlip[facei];
            T u = flip ? fl : fu;
            T l = flip ? fu : fl;

            atomicAdd(cu, u);
            atomicAdd(cl, l);
        }
    }

    __device__
    inline void invoke(float* cD, float* cU, float* cL, const float* fU, const float* fL, const label& facei)
    {
        impl(cD, cU, cL, fU, fL, facei);
    }

    __device__
    inline void invoke(double* cD, double* cU, double* cL, const double* fU, const double* fL, const label& facei)
    {
        #if __CUDA_ARCH__ >= 600
        impl(cD, cU, cL, fU, fL, facei);
        #endif
    }

    __device__
    inline void operator()(const label& fineFacei)
    {
        invoke(coarseDiag, coarseUpper, coarseLower, fineUpper, fineLower, fineFacei);
    }
};


}
}

