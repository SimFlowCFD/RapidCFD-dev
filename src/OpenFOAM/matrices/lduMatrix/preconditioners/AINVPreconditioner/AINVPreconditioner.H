/*---------------------------------------------------------------------------*\
    Sparse Approximate Inverse Preconditioners

    Reference:
    "Algorithm for Sparse Approximate Inverse Preconditioners
    in the Conjugate Gradient Method"

    Ilya B. Labutin, Irina V. Surodina

\*---------------------------------------------------------------------------*/

#ifndef AINVPreconditioner_H
#define AINVPreconditioner_H

#include "lduMatrix.H"

namespace Foam
{

class AINVPreconditioner
:
    public lduMatrix::preconditioner
{
        scalargpuField rD;

        AINVPreconditioner(const AINVPreconditioner&);

        void operator=(const AINVPreconditioner&);

        template<bool normalMult>
        void preconditionImpl
        (
            scalargpuField& w,
            const scalargpuField& r,
            const direction cmpt
        ) const;

        template<bool normalMult, bool useTexture>
        void preconditionImpl
        (
            scalargpuField& w,
            const scalargpuField& r,
            const direction cmpt
        ) const;

public:

    TypeName("AINV");

    AINVPreconditioner
    (
        const lduMatrix::solver&,
        const dictionary& solverControlsUnused
    );

    virtual ~AINVPreconditioner()
    {}

        virtual void precondition
        (
            scalargpuField& wA,
            const scalargpuField& rA,
            const direction cmpt=0
        ) const
        {
            preconditionImpl<true>(wA, rA, cmpt);
        }

        virtual void preconditionT
        (
            scalargpuField& wT,
            const scalargpuField& rT,
            const direction cmpt=0
        ) const
        {
            preconditionImpl<false>(wT, rT, cmpt);
        }
};


} 

#endif
