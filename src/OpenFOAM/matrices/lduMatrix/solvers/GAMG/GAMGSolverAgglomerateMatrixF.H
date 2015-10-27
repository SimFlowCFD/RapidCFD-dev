#pragma once

#include "gpuAtomic.H"

namespace Foam
{

struct luGAMGNonNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x>=0;
    }
};

struct luGAMGNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x<0;
    }
};

struct GAMGSolverAgglomerateAsymmetricFunctor
{
    scalar* coarseUpper;
    scalar* coarseLower;
    scalar* coarseDiag;

    GAMGSolverAgglomerateAsymmetricFunctor
    (
        scalar* _coarseUpper,
        scalar* _coarseLower,
        scalar* _coarseDiag
    ):
        coarseUpper(_coarseUpper),
        coarseLower(_coarseLower),
        coarseDiag(_coarseDiag)
    {}
    
    template<class Tuple>
    __device__
    void operator()(const Tuple& t)
    {
        const label& cFace = thrust::get<0>(t);
        const bool& faceFlipMap = thrust::get<1>(t);
        const scalar& fineUpper = thrust::get<2>(t);
        const scalar& fineLower = thrust::get<3>(t);

        if (cFace >= 0)
        {
            if (!faceFlipMap)
            {
                atomicAdd(coarseUpper+cFace,fineUpper);
                atomicAdd(coarseLower+cFace,fineLower);
            }
            else
            {
                atomicAdd(coarseUpper+cFace,fineLower);
                atomicAdd(coarseLower+cFace,fineUpper);
            }
        }
        else
        {
            atomicAdd(coarseDiag+(-1 - cFace),fineUpper+fineLower);
        }
    }
};

struct GAMGSolverAgglomerateDiagSymmetricFunctor
{
    scalar* coarseDiag;

    GAMGSolverAgglomerateDiagSymmetricFunctor
    (
        scalar* _coarseDiag
    ):
        coarseDiag(_coarseDiag)
    {}
    
    template<class Tuple>
    __device__
    void operator()(const Tuple& t)
    {
        const label& cFace = thrust::get<0>(t);

        if (cFace < 0)
        {
            const scalar& fineUpper = thrust::get<1>(t);
            atomicAdd(coarseDiag+(-1 - cFace),2*fineUpper);
        }
    }
};

struct GAMGSolverAgglomerateSymmetricFunctor
{
    scalar* coarseUpper;

    GAMGSolverAgglomerateSymmetricFunctor
    (
        scalar* _coarseUpper
    ):
        coarseUpper(_coarseUpper)
    {}
    
    template<class Tuple>
    __device__
    void operator()(const Tuple& t)
    {
        const label& cFace = thrust::get<0>(t);

        if (cFace >= 0)
        {
            const scalar& fineUpper = thrust::get<1>(t);
            atomicAdd(coarseUpper+cFace,fineUpper);
        }
    }
};

}

