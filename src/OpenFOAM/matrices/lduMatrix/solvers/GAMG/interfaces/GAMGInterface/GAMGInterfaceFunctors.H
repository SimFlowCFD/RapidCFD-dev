#pragma once

#include "GAMGInterface.H"

namespace Foam
{

template<class Fun, class EqOp>
struct GAMGInterfaceFunctor
{
    const label* cellStart;
    const label* cellFaces;
    const label* interfaceCells;
    Fun fun;
    EqOp op;

    GAMGInterfaceFunctor
    (
        const label* _cellStart,
        const label* _cellFaces,
        const label* _interfaceCells,
        const Fun _fun,
        const EqOp _op
    ):
        cellStart(_cellStart),
        cellFaces(_cellFaces),
        interfaceCells(_interfaceCells),
        fun(_fun),
        op(_op)
    {}

    __device__
    scalar operator()(const label& id,const scalar& s)
    {
        scalar out = s;

        label start = cellStart[id];
        label end = cellStart[id+1];

        for(label i = start; i < end; i++)
        {
            label face = cellFaces[i];
            op(out,fun(interfaceCells[id],face));
        }

        return out;
    }
};

struct GAMGUpdateInterfaceMatrixFunctor
{
    const scalar* coeffs;
    const scalar* pnf;

    GAMGUpdateInterfaceMatrixFunctor
    (
        const scalar* _coeffs,
        const scalar* _pnf
    ):
        coeffs(_coeffs),
        pnf(_pnf)
    {}

    __device__
    scalar operator()(const label& cell, const label& face)
    {
        return coeffs[face]*pnf[face];
    }
};

template<class Fun, class EqOp>
void GAMGInterfaceOperation
(
    gpuList<scalar>& out,
    const GAMGInterface& inter,
    Fun f,
    EqOp o
)
{
    const labelgpuList& cellsSort = inter.sortCells();

    const labelgpuList& cellFaces = inter.cellFaces();
    const labelgpuList& cellFacesStart = inter.cellFacesStart();

    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+cellsSort.size(),
        thrust::make_permutation_iterator
        (
            out.begin(),
            cellsSort.begin()
        ),
        thrust::make_permutation_iterator
        (
            out.begin(),
            cellsSort.begin()
        ),
        GAMGInterfaceFunctor<Fun,EqOp>
        (
            cellFacesStart.data(),
            cellFaces.data(),
            cellsSort.data(),
            f,
            o
        )
    );
}

template<class EqOp>
void GAMGInterfaceOperation
(
    scalargpuList& out,
    const scalargpuList& coeffs,
    const scalargpuList& pnf,
    const GAMGInterface& inter,
    EqOp op
)
{
    GAMGInterfaceOperation
    (
        out,
        inter,
        GAMGUpdateInterfaceMatrixFunctor
        (
            coeffs.data(),
            pnf.data()
        ),
        op
    );
}

inline void GAMGUpdateInterfaceMatrix
(
    scalargpuList& out,
    const scalargpuList& coeffs,
    const scalargpuList& pnf,
    const GAMGInterface& inter,
    const bool reverse
)
{
    if(reverse)
        GAMGInterfaceOperation(out, coeffs, pnf, inter, plusEqOp<scalar>());
    else
        GAMGInterfaceOperation(out, coeffs, pnf, inter, minusEqOp<scalar>());
}

}
