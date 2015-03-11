#pragma once

namespace Foam
{

struct sumlPhiMULESFunctor
{
    const label* own;
    const label* nei;
    const label* ownStart;
    const label* neiStart;
    const label* losort;

    const scalar* lambdaIf;
    const scalar* phiCorrIf;

    scalar* sumlPhip;
    scalar* mSumlPhim;

    sumlPhiMULESFunctor
    (
        const label* _own,
        const label* _nei,
        const label* _ownStart,
        const label* _neiStart,
        const label* _losort,

        const scalar* _lambdaIf,
        const scalar* _phiCorrIf,

        scalar* _sumlPhip,
        scalar* _mSumlPhim
    ):
        own(_own),
        nei(_nei),
        ownStart(_ownStart),
        neiStart(_neiStart),
        losort(_losort),

        lambdaIf(_lambdaIf),
        phiCorrIf(_phiCorrIf),

        sumlPhip(_sumlPhip),
        mSumlPhim(_mSumlPhim)
    {}

    __HOST____DEVICE__
    void operator()(const label& id)
    {
        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;
		
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        scalar sumlPhipTmp = sumlPhip[id];
        scalar mSumlPhimTmp = mSumlPhim[id];

        for(label i = 0; i<oSize; i++)
        {
            label face = oStart + i;

            scalar lambdaPhiCorrf = lambdaIf[face]*phiCorrIf[face];

            if (lambdaPhiCorrf > 0.0)
            {
                sumlPhipTmp += lambdaPhiCorrf;
            }
            else
            {
                mSumlPhimTmp -= lambdaPhiCorrf;
            }
        }

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];

            scalar lambdaPhiCorrf = lambdaIf[face]*phiCorrIf[face];

            if (lambdaPhiCorrf > 0.0)
            {
                mSumlPhimTmp += lambdaPhiCorrf;
            }
            else
            {
                sumlPhipTmp -= lambdaPhiCorrf;
            }
        }

        sumlPhip[id]  = sumlPhipTmp;
        mSumlPhim[id] = mSumlPhimTmp;
    }
};


struct patchSumlPhiMULESFunctor
{
    const label* neiStart;
    const label* losort;
    const label* pCell;

    const scalar* lambdaPf;
    const scalar* phiCorrfPf;

    scalar* sumlPhip;
    scalar* mSumlPhim;

    patchSumlPhiMULESFunctor
    (
        const label* _neiStart,
        const label* _losort,
        const label* _pCell,

        const scalar* _lambdaPf,
        const scalar* _phiCorrfPf,

        scalar* _sumlPhip,
        scalar* _mSumlPhim
    ):
        neiStart(_neiStart),
        losort(_losort),
        pCell(_pCell),

        lambdaPf(_lambdaPf),
        phiCorrfPf(_phiCorrfPf),

        sumlPhip(_sumlPhip),
        mSumlPhim(_mSumlPhim)
    {}

    __HOST____DEVICE__
    void operator()(const label& id)
    {
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        label globalId = pCell[id];

        scalar sumlPhipTmp = sumlPhip[globalId];
        scalar mSumlPhimTmp = mSumlPhim[globalId];

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];

            scalar lambdaPhiCorrf = lambdaPf[face]*phiCorrfPf[face];

            if (lambdaPhiCorrf > 0.0)
            {
                sumlPhipTmp += lambdaPhiCorrf;
            }
            else
            {
                mSumlPhimTmp -= lambdaPhiCorrf;
            }
        }

        sumlPhip[globalId] = sumlPhipTmp;
        mSumlPhim[globalId] = mSumlPhimTmp;
    }
};


template<bool sum>
struct sumlPhipFinalMULESFunctor
{
    __HOST____DEVICE__
    scalar operator()(const scalar& sumlPhi, const thrust::tuple<scalar,scalar>& t)
    {
        scalar psiMinMax = thrust::get<0>(t);
        scalar sumPhi = thrust::get<1>(t);

       return max(min
              (
                    (sumlPhi + psiMinMax)
                   /(sumPhi + (sum?SMALL:-SMALL)),
                    1.0), 0.0
              );
    }
};

struct lambdaIfMULESFunctor
{
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& phiCorrIf,const Tuple& t)
    {
        const scalar& lambdaIf = thrust::get<0>(t);

        if (phiCorrIf > 0.0)
        {
            return min
                (
                    lambdaIf,
                    min(thrust::get<1>(t), thrust::get<2>(t))
                );
        }
        else
        {
            return min
            (
                lambdaIf,
                min(thrust::get<3>(t), thrust::get<4>(t))
            );
        }
    }

};


struct coupledPatchLambdaPfMULESFunctor
{
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& phiCorrfPf, const Tuple& t)
    {
        const scalar& lambdaPf = thrust::get<0>(t);
        
        if (phiCorrfPf > 0.0)
        {
            return min(lambdaPf, thrust::get<1>(t));
        }
        else
        {
            return min(lambdaPf, thrust::get<2>(t));
        }
    }
};

}
