#pragma once

namespace Foam
{
    template<bool fast,int nUnroll>
    struct AINVPreconditionerFunctor 
    {
        const textures<scalar> psi;
        const textures<scalar> rD;
        const scalar* lower;
        const scalar* upper;
        const label* own;
        const label* nei;
        const label* ownStart;
        const label* losortStart;
        const label* losort;

        AINVPreconditionerFunctor
        (
            const textures<scalar> _psi, 
            const textures<scalar> _rD, 
            const scalar* _lower,
            const scalar* _upper,
            const label* _own,
            const label* _nei,
            const label* _ownStart,
            const label* _losortStart,
            const label* _losort
        ):
            psi(_psi),
            rD(_rD),
            lower(_lower),
            upper(_upper),
            own(_own),
            nei(_nei),
            ownStart(_ownStart),
            losortStart(_losortStart),
            losort(_losort)
        {}

        __device__
        scalar operator()(const label& id)
        {
            scalar out = 0;
            scalar tmpSum[2*nUnroll] = {};
            
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;
            
            label nStart = losortStart[id];
            label nSize = losortStart[id+1] - nStart;

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<oSize)
                {
                    label face = oStart + i;

                    tmpSum[i] = upper[face]*rD[nei[face]]*psi[nei[face]];
                }
            }

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<nSize)
                {
                    label face = nStart + i;
                    if(!fast)
                         face = losort[face];

                    tmpSum[i+nUnroll] = lower[face]*rD[own[face]]*psi[own[face]];
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

                out += upper[face]*rD[nei[face]]*psi[nei[face]];
            }
            
            #pragma unroll 2
            for(label i = nUnroll; i<nSize; i++)
            {
                 label face = nStart + i;
                 if(!fast)
                      face = losort[face];

                 out += lower[face]*rD[own[face]]*psi[own[face]];
            }

            
            return rD[id]*(psi[id]-out);
        }
    };
}
