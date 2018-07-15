#pragma once

namespace Foam
{

    template<bool fast,int nUnroll>
    struct JacobiSmootherFunctor 
    {
        const scalar* psi;
        const scalar* diag;
        const scalar* b;
        const scalar* lower;
        const scalar* upper;
        const __restrict__ label* own;
        const __restrict__ label* nei;
        const __restrict__ label* ownStart;
        const __restrict__ label* losortStart;
        const __restrict__ label* losort;
        const __restrict__ scalar omega;

        JacobiSmootherFunctor
        (
            scalar _omega,
            const scalar* _psi,
            const scalar* _diag, 
            const scalar* _b, 
            const scalar* _lower,
            const scalar* _upper,
            const __restrict__ label* _own,
            const __restrict__ label* _nei,
            const __restrict__ label* _ownStart,
            const __restrict__ label* _losortStart,
            const __restrict__ label* _losort
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

        #if __CUDA_ARCH__ >= 350
        
            scalar out = 0;
            scalar tmpSum[2*nUnroll] = {};
            const scalar rD = 1.0/__ldg(&diag[id]);

            scalar extra = (1 - omega)*__ldg(&psi[id]) + omega*rD*__ldg(&b[id]);
            
            label oStart = __ldg(&ownStart[id]);
            label oSize = __ldg(&ownStart[id+1]) - oStart;
            
            label nStart = __ldg(&losortStart[id]);
            label nSize = __ldg(&losortStart[id+1]) - nStart;

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<oSize)
                {
                    label face = oStart + i;
                    tmpSum[i] = __ldg(&upper[face])*__ldg(&psi[nei[face]]);
                }
            }

            for(label i = 0; i<nUnroll; i++)
            {
                if(i<nSize)
                {
                     label face = nStart + i;
                     if( ! fast) 
                         face = __ldg(&losort[face]);

                     tmpSum[i+nUnroll] = __ldg(&lower[face])*__ldg(&psi[own[face]]); 
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
                out += __ldg(&upper[face])*__ldg(&psi[nei[face]]);
            }
            
            #pragma unroll 2
            for(label i = nUnroll; i<nSize; i++)
            {
                 label face = nStart + i;
                 if( ! fast) 
                     face = __ldg(&losort[face]);

                 out += __ldg(&lower[face])*__ldg(&psi[own[face]]);
            }
            
            return extra - omega*rD*out;
 
        #else
    
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
 
        #endif
 
        }
    };

}
