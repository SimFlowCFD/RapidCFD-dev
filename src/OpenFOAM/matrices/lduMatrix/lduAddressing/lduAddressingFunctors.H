#ifndef lduAddressingFunctors_H
#define lduAddressingFunctors_H

#include "lduAddressing.H"
#include "ops.H"

namespace Foam
{

template<class Type, class OwnFun, class NeiFun, class OwnOp, class NeiOp>
struct lduAddressingFunctor
{
    const label* ownStart;
    const label* neiStart;
    const label* losort;
    OwnFun oFun;
    NeiFun nFun;
    OwnOp oOp;
    NeiOp nOp;

    lduAddressingFunctor
    (
        const label* _ownStart,
        const label* _neiStart,
        const label* _losort,
        const OwnFun _oFun,
        const NeiFun _nFun,
        const OwnOp _oOp,
        const NeiOp _nOp
    ):
        ownStart(_ownStart),
        neiStart(_neiStart),
        losort(_losort),
        oFun(_oFun),
        nFun(_nFun),
        oOp(_oOp),
        nOp(_nOp)
    {}

    __HOST____DEVICE__
    Type operator()(const label& id,const Type& d)
    {
        Type out = d;
        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;
		
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<oSize; i++)
        {
            label face = oStart + i;
            out = oOp(out,oFun(id,face));
        }

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            out = nOp(out,nFun(id,face));
        }

        return out;
    }
};

#define MAX_NEI_SIZE 3
template<class OwnFun, class NeiFun>
struct lduAddressingFastFunctor
{
    const label* ownStart;
    const label* neiStart;
    OwnFun oFun;
    NeiFun nFun;

    lduAddressingFastFunctor
    (
        const label* _ownStart,
        const label* _neiStart,
        const OwnFun _oFun,
        const NeiFun _nFun
    ):
        ownStart(_ownStart),
        neiStart(_neiStart),
        oFun(_oFun),
        nFun(_nFun)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& id,const scalar& d)
    {
        scalar out = d;
        scalar tmpSum[2*MAX_NEI_SIZE] = {};
        scalar nExtra = 0;

        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;
		
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<MAX_NEI_SIZE; i++)
        {
            if(i<oSize)
            {
                label face = oStart + i;
                tmpSum[i] = oFun(id,face);
            }
        }

        for(label i = 0; i<MAX_NEI_SIZE; i++)
        {
            if(i<nSize)
            {
                 label face = nStart + i;
                   
                 tmpSum[i+MAX_NEI_SIZE] = nFun(id,face);
            }
        }

        for(label i = 0; i<2*MAX_NEI_SIZE; i++)
        {
            out+= tmpSum[i]; 
        }
           
        for(label i = MAX_NEI_SIZE; i<oSize; i++)
        {
            label face = oStart + i;
                
            out += oFun(id,face);
        }
            
        for(label i = MAX_NEI_SIZE; i<nSize; i++)
        {
            label face = nStart + i;

            nExtra += nFun(id,face); 
        }  

        return out;
    }
};
#undef MAX_NEI_SIZE

template<class Type,class Fun,class Op>
struct lduAddressingPatchFunctor
{
    const label* neiStart;
    const label* losort;
    const label* pcells;
    Fun fun;
    Op op;

    lduAddressingPatchFunctor
    (
        const label* _neiStart,
        const label* _losort,
        const label* _pcells,
        const Fun _fun,
        const Op _op
    ):
        neiStart(_neiStart),
        losort(_losort),
        pcells(_pcells),
        fun(_fun),
        op(_op)
    {}

    __HOST____DEVICE__
    Type operator()(const label& id,const Type& s)
    {
        Type out = s;

        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;
        label cellI = pcells[id];

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            out = op(out,fun(cellI,face));
        }

        return out;
    }
};

template<class Type,class LUType,class Op>
struct matrixCoeffsMultiplyFunctor
{
    const Type* psi;
    const LUType* coeffs;
    const label* addr;
    Op op;

    matrixCoeffsMultiplyFunctor
    (
        const Type* _psi,
        const LUType* _coeffs,
        const label* _addr,
        const Op _op
    ):
        psi(_psi),
        coeffs(_coeffs),
        addr(_addr),
        op(_op)
    {}

    __HOST____DEVICE__
    Type operator()(const label& cell, const label& face) 
    {
        return op(coeffs[face]*psi[addr[face]]);
    }
};

template<class Type,class Op>
struct matrixCoeffsFunctor
{
    const Type* coeffs;
    Op op;

    matrixCoeffsFunctor
    (
        const Type* _coeffs,
        const Op _op
    ):
        coeffs(_coeffs),
        op(_op)
    {}

    __HOST____DEVICE__
    Type operator()(const label& cell, const label& face)
    {
        return op(coeffs[face]);
    }
};

template<class Type>
struct matrixInterfaceFunctor
{
    const scalar* coeffs;
    const Type* val;

    matrixInterfaceFunctor
    (
        const scalar* _coeffs, 
        const Type* _val
    ):
        coeffs(_coeffs),
        val(_val)
    {}

    __HOST____DEVICE__
    Type operator()(const label&, const label& id)
    {
        return -coeffs[id]*val[id];
    }
};

template<class Type,class Input, class OwnFun,class NeiFun,class OwnOp,class NeiOp>
inline void matrixOperation
(
    Input in, 
    gpuList<Type>& out,
    const lduAddressing& addr,
    OwnFun o, 
    NeiFun n,
    OwnOp oo,
    NeiOp no
)
{
    const labelgpuList& ownStart = addr.ownerStartAddr();
    const labelgpuList& losortStart = addr.losortStartAddr();
    const labelgpuList& losort = addr.losortAddr();
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+addr.size(),
        in,
        out.begin(),
        lduAddressingFunctor<Type,OwnFun,NeiFun,OwnOp,NeiOp>
        (
            ownStart.data(),
            losortStart.data(),
            losort.data(),
            o,
            n,
            oo,
            no
        )
    );
}


template<class Type,class Input, class OwnFun,class NeiFun>
inline void matrixOperation
(
    Input in, 
    gpuList<Type>& out,
    const lduAddressing& addr,
    OwnFun o, 
    NeiFun n
)
{
    matrixOperation
    (
        in, 
        out,
        addr,
        o, 
        n,
        sumOp<Type>(),
        sumOp<Type>()
    );
}

template<class Input, class OwnFun,class NeiFun>
inline void matrixFastOperation
(
    Input in, 
    gpuList<scalar>& out,
    const lduAddressing& addr,
    OwnFun o, 
    NeiFun n
)
{
    const labelgpuList& ownStart = addr.ownerStartAddr();
    const labelgpuList& losortStart = addr.losortStartAddr();
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+addr.size(),
        in,
        out.begin(),
        lduAddressingFastFunctor<OwnFun,NeiFun>
        (
            ownStart.data(),
            losortStart.data(),
            o,
            n
        )
    );
}

template<class Type, class Fun, class Op>
inline void matrixPatchOperation
(
    label patchI, 
    gpuList<Type>& out,
    const lduAddressing& addr,
    Fun f,
    Op o
)
{
    const labelgpuList& pcells = addr.patchSortCells(patchI);

    const labelgpuList& losort = addr.patchSortAddr(patchI);
    const labelgpuList& losortStart = addr.patchSortStartAddr(patchI);
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+pcells.size(),
        thrust::make_permutation_iterator(out.begin(),pcells.begin()),
        thrust::make_permutation_iterator(out.begin(),pcells.begin()),
        lduAddressingPatchFunctor<Type,Fun,Op>
        (
            losortStart.data(),
            losort.data(),
            pcells.data(),
            f,
            o
        )
    );
}

template<class Type, class Fun>
inline void matrixPatchOperation
(
    label patchI, 
    gpuList<Type>& out,
    const lduAddressing& addr,
    Fun f
)
{
    matrixPatchOperation
    (
        patchI, 
        out,
        addr,
        f,
        sumOp<Type>()
    );
}

}

#endif
