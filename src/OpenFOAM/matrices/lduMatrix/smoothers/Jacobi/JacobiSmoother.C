#include "JacobiSmoother.H"
#include "textureConfig.H"

namespace Foam
{
    defineTypeNameAndDebug(JacobiSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<JacobiSmoother>
        addJacobiSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<JacobiSmoother>
        addJacobiSmootherAsymMatrixConstructorToTable_;   

    #define MAX_NEI_SIZE 3
	
    template<bool useTexture>
    struct JacobiSmootherFunctor 
    {
        const scalar* psi;
        const scalar* diag;
        const scalar* b;
        const scalar* lower;
        const scalar* upper;
        const label* own;
        const label* nei;
        const label* losort;
        const label* ownStart;
        const label* losortStart;
        const scalar omega;

        JacobiSmootherFunctor
        (
            scalar _omega,
            const scalar* _psi, 
            const scalar* _diag, 
            const scalar* _b, 
            const scalar* _lower,
            const scalar* _upper,
            const label* _own,
            const label* _nei,
            const label* _losort,
            const label* _ownStart,
            const label* _losortStart
        ):
             psi(_psi),
             diag(_diag),
             b(_b),
             lower(_lower),
             upper(_upper),
             own(_own),
             nei(_nei),
             losort(_losort),
             ownStart(_ownStart),
             losortStart(_losortStart),
             omega(_omega)
        {}

        __HOST____DEVICE__
        scalar operator()(const label& id)
        {
            scalar out = 0;
            scalar tmpSum[2*MAX_NEI_SIZE] = {};
            const scalar rD = 1.0/diag[id];

            scalar extra = (1 - omega)*psi[id] + omega*rD*b[id];
            
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;
            
            label nStart = losortStart[id];
            label nSize = losortStart[id+1] - nStart;

            for(label i = 0; i<MAX_NEI_SIZE; i++)
            {
                if(i<oSize)
                {
                    label face = oStart + i;
                    tmpSum[i] = upper[face]*fetch<useTexture>(nei[face],psi);
                }
            }

            for(label i = 0; i<MAX_NEI_SIZE; i++)
            {
                if(i<nSize)
                {
                     label face = losort[nStart + i];
                     tmpSum[i+MAX_NEI_SIZE] = lower[face]*fetch<useTexture>(own[face],psi); 
                }
            }

            for(label i = 0; i<2*MAX_NEI_SIZE; i++)
            {
                out+= tmpSum[i]; 
            }
            
            for(label i = MAX_NEI_SIZE; i<oSize; i++)
            {
                label face = oStart + i;
                out += upper[face]*fetch<useTexture>(nei[face],psi);
            }
            
            
            for(label i = MAX_NEI_SIZE; i<nSize; i++)
            {
                 label face = losort[nStart + i];

                 out += lower[face]*fetch<useTexture>(own[face],psi);
            }

            
            return extra - omega*rD*out;
        }
    };
    
    #undef MAX_NEI_SIZE
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
    scalargpuField Apsi(psi.size());
    scalargpuField sourceTmp(source.size());


    const labelgpuList& l = matrix_.lduAddr().lowerAddr();
    const labelgpuList& u = matrix_.lduAddr().upperAddr();
    const labelgpuList& losort = matrix_.lduAddr().losortAddr();

    const labelgpuList& ownStart = matrix_.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = matrix_.lduAddr().losortStartAddr();

    const scalargpuField& Lower = matrix_.lower();
    const scalargpuField& Upper = matrix_.upper();
    const scalargpuField& Diag = matrix_.diag();

    const bool textureCanBeUsed = psi.size() > TEXTURE_MINIMUM_SIZE;


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
        sourceTmp = source;

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

        if(textureCanBeUsed)
        {

        bind(psi.data());

        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+psi.size(),
            Apsi.begin(),
            JacobiSmootherFunctor<true>
            (
                omega_,
                psi.data(),
                Diag.data(),
                sourceTmp.data(),
                Lower.data(),
                Upper.data(),
                l.data(),
                u.data(),
                losort.data(),
                ownStart.data(),
                losortStart.data()
            )
        );

        unbind(psi.data());

        }
        else
        {

        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+psi.size(),
            Apsi.begin(),
            JacobiSmootherFunctor<false>
            (
                omega_,
                psi.data(),
                Diag.data(),
                sourceTmp.data(),
                Lower.data(),
                Upper.data(),
                l.data(),
                u.data(),
                losort.data(),
                ownStart.data(),
                losortStart.data()
            )
        );

        }

        psi = Apsi;
    }

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }
}

