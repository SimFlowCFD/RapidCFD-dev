#ifndef JacobiSmoother_H
#define JacobiSmoother_H

#include "lduMatrix.H"

namespace Foam
{

class JacobiSmoother
:
    public lduMatrix::smoother
{
    scalar omega_;

public:

    TypeName("Jacobi");

    JacobiSmoother
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<gpuField, scalar>& interfaceBouCoeffs,
        const FieldField<gpuField, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual void smooth
    (
        scalargpuField& psi,
        const scalargpuField& source,
        const direction cmpt,
        const label nSweeps
    ) const;
};

}

#endif

