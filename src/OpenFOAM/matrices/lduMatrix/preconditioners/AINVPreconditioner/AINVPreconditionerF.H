#pragma once

namespace Foam
{
    template<int nUnroll>
    struct AINVPreconditionerFunctor 
    {
        textures<scalar> psi;
        textures<scalar> rD;
        const scalar* lower;
        const scalar* upper;
        const label* own;
        const label* nei;
        const label* ownStart;
        const label* losortStart;

        AINVPreconditionerFunctor
        (
            textures<scalar> _psi, 
            textures<scalar> _rD, 
            const scalar* _lower,
            const scalar* _upper,
            const label* _own,
            const label* _nei,
            const label* _ownStart,
            const label* _losortStart
        ):
            psi(_psi),
            rD(_rD),
            lower(_lower),
            upper(_upper),
            own(_own),
            nei(_nei),
            ownStart(_ownStart),
            losortStart(_losortStart)
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

                    tmpSum[i+nUnroll] = lower[face]*rD[own[face]]*psi[own[face]];
                }
            }

            for(label i = 0; i<2*nUnroll; i++)
            {
                out+= tmpSum[i]; 
            }
            
            for(label i = nUnroll; i<oSize; i++)
            {
                label face = oStart + i;

                out += upper[face]*rD[nei[face]]*psi[nei[face]];
            }
            
            
            for(label i = nUnroll; i<nSize; i++)
            {
                 label face = nStart + i;

                 out += lower[face]*rD[own[face]]*psi[own[face]];
            }

            
            return rD[id]*(psi[id]-out);
        }
    };
}
