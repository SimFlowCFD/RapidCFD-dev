#ifndef lduMatrixFunctors_H
#define lduMatrixFunctors_H

namespace Foam
{
    struct lduMatrixDiagonalFunctor : public std::unary_function<thrust::tuple<scalar,scalar>,scalar>{
        __host__ __device__
        scalar operator()(const thrust::tuple<scalar,scalar>& c){
            return thrust::get<0>(c) * thrust::get<1>(c); 
        }
    };
    
    struct lduMatrixDiagonalResidualFunctor : public std::unary_function<thrust::tuple<scalar,scalar,scalar>,scalar>{
        __host__ __device__
        scalar operator()(const thrust::tuple<scalar,scalar,scalar>& c){
            return thrust::get<0>(c) -          //source 
                           thrust::get<1>(c) *  //diagonal
                           thrust::get<2>(c);   //flux
        }
    };

    template<class Type,bool addOffDiagonal, bool normalMult>
    struct lduMatrixMultiplyFunctor : public std::binary_function<scalar,label,scalar>{
        const Type* psi;
        const scalar* lower;
        const scalar* upper;
        const label* own;
        const label* nei;
        const label* losort;
        lduMatrixMultiplyFunctor(const Type* _psi, 
                                 const scalar* _lower,
                                 const scalar* _upper,
                                 const label* _own,
                                 const label* _nei,
                                 const label* _losort):
             psi(_psi),
             lower(_lower),
             upper(_upper),
             own(_own),
             nei(_nei),
             losort(_losort)
        {}
        __HOST____DEVICE__
        Type operator()(const Type& d,const thrust::tuple<label,label,label,label>& t){
            Type out = d;
            label oStart = thrust::get<0>(t);
            label oSize = thrust::get<1>(t) - oStart;
            
            label nStart = thrust::get<2>(t);
            label nSize = thrust::get<3>(t) - nStart;

            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                if(addOffDiagonal)
                {
                    if(normalMult)
                        out += upper[face]*psi[nei[face]]; 
                    else
                        out += lower[face]*psi[nei[face]]; 
                }
                else
                {
                    if(normalMult)
                        out -= upper[face]*psi[nei[face]]; 
                    else
                        out -= lower[face]*psi[nei[face]]; 
                }
            }


            for(label i = 0; i<nSize; i++)
            {
                 label face = losort[nStart + i];
                if(addOffDiagonal)
                {
                    if(normalMult)
                        out += lower[face]*psi[own[face]]; 
                    else
                        out += upper[face]*psi[own[face]];
                }
                else
                {
                    if(normalMult)
                        out -= lower[face]*psi[own[face]]; 
                    else
                        out -= upper[face]*psi[own[face]]; 
                }
            }

            return out;
        }
    };

    template<bool addOffDiagonal>
    struct lduMatrixSumFunctor : public std::binary_function<scalar,label,scalar>{
        const scalar* lower;
        const scalar* upper;
        const label* ownStart;
        const label* neiStart;
        const label* own;
        const label* nei;
        const label* losort;
        lduMatrixSumFunctor(const scalar* _lower,
                                 const scalar* _upper,
                                 const label* _ownStart,
                                 const label* _neiStart,
                                 const label* _own,
                                 const label* _nei,
                                 const label* _losort):
             lower(_lower),
             upper(_upper),
             ownStart(_ownStart),
             neiStart(_neiStart),
             own(_own),
             nei(_nei),
             losort(_losort)
        {}
        __HOST____DEVICE__
        scalar operator()(const scalar& diag,const label& id){
            scalar out = diag;
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;

            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                if(addOffDiagonal)
                {
                    out += upper[face];
                }
                else
                {
                    out -= upper[face];
                }
            }

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                if(addOffDiagonal)
                {
                    out += lower[face];
                }
                else
                {
                    out -= lower[face]; 
                }
            }

            return out;
        }
    };

    struct lduMatrixPatchSubtractFunctor : public std::binary_function<label,scalar,scalar>{
        const scalar* coeff;
        const label* neiStart;
        const label* losort;
        lduMatrixPatchSubtractFunctor(const scalar* _coeff,
                                 const label* _neiStart,
                                 const label* _losort):
             coeff(_coeff),
             neiStart(_neiStart),
             losort(_losort)
        {}
        __HOST____DEVICE__
        scalar operator()(const label& id,const scalar& s){
            scalar out = s;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out -= coeff[face];
            }

            return out;
        }
    };

    struct psiPlusAlphaYAFunctor
    {
        const scalar alpha;

        psiPlusAlphaYAFunctor(scalar _alpha): alpha(_alpha) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& psi, const scalar& yA)
        {
                return psi + alpha*yA;
        }
    };

    struct rAPlusAlphaPAFunctor
    {
        const scalar alpha;

        rAPlusAlphaPAFunctor(scalar _alpha): alpha(_alpha) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& rA, const scalar& pA)
        {
                return rA + alpha*pA;
        }
    };

    struct rAMinusAlphaPAFunctor
    {
        const scalar alpha;

        rAMinusAlphaPAFunctor(scalar _alpha): alpha(_alpha) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& rA, const scalar& pA)
        {
                return rA - alpha*pA;
        }
    };

    struct rAPlusBetaPAFunctor
    {
        const scalar beta;

        rAPlusBetaPAFunctor(scalar _beta): beta(_beta) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& rA, const scalar& pA)
        {
                return rA + beta*pA;
        }
    };

    struct rAMinusAlphaAYAFunctor
    {
        const scalar alpha;

        rAMinusAlphaAYAFunctor(scalar _alpha): alpha(_alpha) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& rA, const scalar& AyA)
        {
                return rA - alpha*AyA;
        }
    };

    struct sAMinusOmegaTAFunctor
    {
        const scalar omega;

        sAMinusOmegaTAFunctor(scalar _omega): omega(_omega) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& sA, const scalar& tA)
        {
                return sA - omega*tA;
        }
    };

    struct pAMinusOmegaAYAFunctor
    {
        const scalar omega;

        pAMinusOmegaAYAFunctor(scalar _omega): omega(_omega) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& pA, const scalar& AyA)
        {
                return pA - omega*AyA;
        }
    };

    struct rAPlusBetaResult1Functor
    {
        const scalar beta;

        rAPlusBetaResult1Functor(scalar _beta): beta(_beta) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& rA, const scalar& result1)
        {
                return rA + beta*result1;
        }
    };

    struct psiPlusOmegaZAFunctor
    {
        const scalar omega;

        psiPlusOmegaZAFunctor(scalar _omega): omega(_omega) {}

        __HOST____DEVICE__
        scalar operator()(const scalar& psi, const scalar& zA)
        {
                return psi + omega*zA;
        }
    };

}

#endif
