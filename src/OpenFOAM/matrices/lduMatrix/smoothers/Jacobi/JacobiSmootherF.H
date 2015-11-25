#pragma once

#include "textures.H"

namespace Foam
{

    template<bool fast,int nUnroll>
    struct JacobiSmootherFunctor 
    {
        const textures<scalar> psi;
        const scalar* diag;
        const scalar* b;
        const scalar* lower;
        const scalar* upper;
        const label* own;
        const label* nei;
        const label* ownStart;
        const label* losortStart;
        const label* losort;
        const scalar omega;

        JacobiSmootherFunctor
        (
            scalar _omega,
            const textures<scalar> _psi, 
            const scalar* _diag, 
            const scalar* _b, 
            const scalar* _lower,
            const scalar* _upper,
            const label* _own,
            const label* _nei,
            const label* _ownStart,
            const label* _losortStart,
            const label* _losort
        ):
            psi(_psi),
            diag(_diag),
            b(_b),
            lower(_lower),
            upper(_upper),
            own(_own),
            nei(_nei),
            ownStart(_ownStart),
            losortStart(_losortStart),
            losort(_losort),
            omega(_omega)
        {}

        __device__
        scalar operator()(const label& id)
        {
            scalar out = 0;
            scalar tmpSum[2*nUnroll] = {};
            const scalar rD = 1.0/diag[id];

            scalar extra = (1 - omega)*psi[id] + omega*rD*b[id];
            
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;
            
            label nStart = losortStart[id];
            label nSize = losortStart[id+1] - nStart;

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<oSize)
                {
                    label face = oStart + i;
                    tmpSum[i] = upper[face]*psi[nei[face]];
                }
            }

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<nSize)
                {
                     label face = nStart + i;
                     if( ! fast) 
                         face = losort[face];

                     tmpSum[i+nUnroll] = lower[face]*psi[own[face]]; 
                }
            }

            #pragma unroll
            for(label i = 0; i<2*nUnroll; i++)
            {
                out+= tmpSum[i]; 
            }
            
            #pragma unroll 2
            for(label i = nUnroll; i<oSize; i++)
            {
                label face = oStart + i;
                out += upper[face]*psi[nei[face]];
            }
            
            #pragma unroll 2
            for(label i = nUnroll; i<nSize; i++)
            {
                 label face = nStart + i;
                 if( ! fast) 
                     face = losort[face];

                 out += lower[face]*psi[own[face]];
            }
            
            return extra - omega*rD*out;
        }
    };

}
