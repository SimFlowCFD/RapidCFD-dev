#ifndef lduAddressingFunctors_H
#define lduAddressingFunctors_H

#include "lduAddressing.H"

namespace Foam
{

template<class Type,class OwnFun,class NeiFun>
struct lduAddressingFunctor : public std::binary_function<label,Type,Type>
{
    const label* own;
    const label* nei;
    const label* ownStart;
    const label* neiStart;
    const label* losort;
    OwnFun oFun;
    NeiFun nFun;
    lduAddressingFunctor
    (
        const label* _own,
        const label* _nei,
        const label* _ownStart,
        const label* _neiStart,
        const label* _losort,
        const OwnFun _oFun,
        const NeiFun _nFun
    ):
        own(_own),
        nei(_nei),
        ownStart(_ownStart),
        neiStart(_neiStart),
        losort(_losort),
        oFun(_oFun),
        nFun(_nFun)
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
            out += oFun(face);
        }

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            out += nFun(face);
        }

        return out;
    }
};

template<class Type,class Fun>
struct lduAddressingPatchFunctor : public std::binary_function<label,Type,Type>
{
    const label* neiStart;
    const label* losort;
    Fun fun;
    lduAddressingPatchFunctor(const label* _neiStart,const label* _losort,const Fun _fun):
                              neiStart(_neiStart),losort(_losort),fun(_fun){}
    __HOST____DEVICE__
    Type operator()(const label& id,const Type& s)
    {
        Type out = s;

        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            out += fun(face);
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
    matrixCoeffsMultiplyFunctor(const Type* _psi,const LUType* _coeffs,const label* _addr,const Op _op):
	                        psi(_psi),coeffs(_coeffs),addr(_addr),op(_op){}
    __HOST____DEVICE__
    Type operator()(const label& face) 
    {
        return op(coeffs[face]*psi[addr[face]]);
    }
};

template<class Type,class Op>
struct matrixCoeffsFunctor
{
    const Type* coeffs;
    Op op;
    matrixCoeffsFunctor(const Type* _coeffs,const Op _op):
	                coeffs(_coeffs),op(_op){}
    __HOST____DEVICE__
    Type operator()(const label& face)
    {
        return op(coeffs[face]);
    }
};


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
    const labelgpuList& l = addr.lowerAddr();
    const labelgpuList& u = addr.upperAddr();
    const labelgpuList& losort = addr.losortAddr();

    const labelgpuList& ownStart = addr.ownerStartAddr();
    const labelgpuList& losortStart = addr.losortStartAddr();
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+addr.size(),
        in,
        out.begin(),
        lduAddressingFunctor<Type,OwnFun,NeiFun>
        (
            l.data(),
            u.data(),
            ownStart.data(),
            losortStart.data(),
            losort.data(),
            o,
            n
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
    const labelgpuList& pcells = addr.patchSortCells(patchI);

    const labelgpuList& losort = addr.patchSortAddr(patchI);
    const labelgpuList& losortStart = addr.patchSortStartAddr(patchI);
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+pcells.size(),
        thrust::make_permutation_iterator(out.begin(),pcells.begin()),
        thrust::make_permutation_iterator(out.begin(),pcells.begin()),
        lduAddressingPatchFunctor<Type,Fun>
        (
            losortStart.data(),
            losort.data(),
            f
        )
    );
}

}

#endif
