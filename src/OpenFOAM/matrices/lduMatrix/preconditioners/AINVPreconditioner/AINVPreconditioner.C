#include "AINVPreconditioner.H"
#include "AINVPreconditionerF.H"
#include "lduMatrixSolutionCache.H"

namespace Foam
{
    defineTypeNameAndDebug(AINVPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<AINVPreconditioner>
        addAINVPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<AINVPreconditioner>
        addAINVPreconditionerAsymMatrixConstructorToTable_;
}

Foam::AINVPreconditioner::AINVPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD
    (
        lduMatrixSolutionCache::first(sol.matrix().diag().size()),
        sol.matrix().diag().size()
    ),
    rDTex(rD)
{ 
  
    const scalargpuField& Diag = solver_.matrix().diag();

    thrust::transform
    (
        Diag.begin(),
        Diag.end(),
        rD.begin(),
        divideOperatorSFFunctor<scalar,scalar,scalar>(1.0)
    );
}

Foam::AINVPreconditioner::~AINVPreconditioner()
{
   rDTex.destroy();
}

template<bool normalMult>
void Foam::AINVPreconditioner::preconditionImpl
(
    scalargpuField& w,
    const scalargpuField& r,
    const direction d
) const
{
    bool fastPath = lduMatrixSolutionCache::favourSpeed;

    const labelgpuList& l = fastPath? 
                            solver_.matrix().lduAddr().ownerSortAddr():
                            solver_.matrix().lduAddr().lowerAddr();
    const labelgpuList& u = solver_.matrix().lduAddr().upperAddr();

    const labelgpuList& ownStart = solver_.matrix().lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = solver_.matrix().lduAddr().losortStartAddr();
    const labelgpuList& losort = solver_.matrix().lduAddr().losortAddr();

    const scalargpuField& Lower = normalMult?
                                  (fastPath?solver_.matrix().lowerSort():solver_.matrix().lower()):
                                  (fastPath?solver_.matrix().upperSort():solver_.matrix().upper());

    const scalargpuField& Upper = normalMult?
                                  solver_.matrix().upper():
                                  solver_.matrix().lower();

    textures<scalar> rTex(r);

    if(fastPath)
    {
        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+r.size(),
            w.begin(),
            AINVPreconditionerFunctor<true,3>
            (
                rTex,
                rDTex,
                Lower.data(),
                Upper.data(),
                l.data(),
                u.data(),
                ownStart.data(),
                losortStart.data(),
                losort.data()
            )
        );
    }
    else
    {
        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+r.size(),
            w.begin(),
            AINVPreconditionerFunctor<false,3>
            (
                rTex,
                rDTex,
                Lower.data(),
                Upper.data(),
                l.data(),
                u.data(),
                ownStart.data(),
                losortStart.data(),
                losort.data()
            )
        );
    }

    rTex.destroy();
}

