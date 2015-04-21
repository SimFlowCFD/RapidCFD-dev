#include "JacobiSmoother.H"
#include "JacobiSmootherF.H"
#include "BasicCache.H"

namespace Foam
{
    defineTypeNameAndDebug(JacobiSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<JacobiSmoother>
        addJacobiSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<JacobiSmoother>
        addJacobiSmootherAsymMatrixConstructorToTable_;   

    class JacobiCache
    {
        static PtrList<scalargpuField> psiCache;
        static PtrList<scalargpuField> sourceCache;

        public:

        static scalargpuField& psi(label level, label size)
        {
            return cache::retrieve(psiCache,level,size);
        }

        static scalargpuField& source(label level, label size)
        {
            return cache::retrieve(sourceCache,level,size);
        }
    };

    PtrList<scalargpuField> JacobiCache::psiCache(1);
    PtrList<scalargpuField> JacobiCache::sourceCache(1);
}

Foam::JacobiSmoother::JacobiSmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const FieldField<gpuField, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    omega_(0.8)
{
    solverControls.readIfPresent("omega", omega_);
}

void Foam::JacobiSmoother::smooth
(
    scalargpuField& psi,
    const scalargpuField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    scalargpuField& Apsi = JacobiCache::psi(matrix_.level(),psi.size());
    scalargpuField& sourceTmp = JacobiCache::source(matrix_.level(),source.size());

    const labelgpuList& l = matrix_.lduAddr().ownerSortAddr();
    const labelgpuList& u = matrix_.lduAddr().upperAddr();

    const labelgpuList& ownStart = matrix_.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = matrix_.lduAddr().losortStartAddr();

    const scalargpuField& Lower = matrix_.lowerSort();
    const scalargpuField& Upper = matrix_.upper();
    const scalargpuField& Diag = matrix_.diag();

    textures<scalar> psiTex(psi);

    FieldField<gpuField, scalar>& mBouCoeffs =
        const_cast<FieldField<gpuField, scalar>&>
        (
            interfaceBouCoeffs_
        );

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        thrust::copy
        (
            source.begin(),
            source.end(),
            sourceTmp.begin()
        );

        matrix_.initMatrixInterfaces
        (
            interfaceBouCoeffs_,
            interfaces_,
            psi,
            sourceTmp,
            cmpt
        );

        matrix_.updateMatrixInterfaces
        (
            interfaceBouCoeffs_,
            interfaces_,
            psi,
            sourceTmp,
            cmpt
        );

        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+psi.size(),
            Apsi.begin(),
            JacobiSmootherFunctor<3>
            (
                omega_,
                psiTex,
                Diag.data(),
                sourceTmp.data(),
                Lower.data(),
                Upper.data(),
                l.data(),
                u.data(),
                ownStart.data(),
                losortStart.data()
            )
        );

        thrust::copy
        (
            Apsi.begin(),
            Apsi.begin()+psi.size(),
            psi.begin()
        );
    }

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }

    psiTex.destroy();
}

