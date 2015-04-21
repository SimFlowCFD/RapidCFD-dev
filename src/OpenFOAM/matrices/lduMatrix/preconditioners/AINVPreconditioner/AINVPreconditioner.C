#include "AINVPreconditioner.H"
#include "AINVPreconditionerF.H"
#include "BasicCache.H"

namespace Foam
{
    defineTypeNameAndDebug(AINVPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<AINVPreconditioner>
        addAINVPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<AINVPreconditioner>
        addAINVPreconditionerAsymMatrixConstructorToTable_;

    class AINVPreconditionerCache
    {
        static PtrList<scalargpuField> rDCache;
        public:

        static scalargpuField& rD(label level, label size)
        {
            return cache::retrieve(rDCache,level,size);
        }
    };

    PtrList<scalargpuField> AINVPreconditionerCache::rDCache(1);
}

Foam::AINVPreconditioner::AINVPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rDPtr
    (
        &AINVPreconditionerCache::rD
        (
            sol.matrix().level(),
            sol.matrix().diag().size()
        )
    ),
    rDTex(*rDPtr)
{ 
  
    const scalargpuField& Diag = solver_.matrix().diag();
    scalargpuField& rD = *rDPtr;

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
    const labelgpuList& l = solver_.matrix().lduAddr().ownerSortAddr();
    const labelgpuList& u = solver_.matrix().lduAddr().upperAddr();

    const labelgpuList& ownStart = solver_.matrix().lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = solver_.matrix().lduAddr().losortStartAddr();

    const scalargpuField& Lower = normalMult?
                                  solver_.matrix().lowerSort():
                                  solver_.matrix().upperSort();

    const scalargpuField& Upper = normalMult?
                                  solver_.matrix().upper():
                                  solver_.matrix().lower();

    textures<scalar> rTex(r);

    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+r.size(),
        w.begin(),
        AINVPreconditionerFunctor<3>
        (
            rTex,
            rDTex,
            Lower.data(),
            Upper.data(),
            l.data(),
            u.data(),
            ownStart.data(),
            losortStart.data()
        )
    );

    rTex.destroy();
}

