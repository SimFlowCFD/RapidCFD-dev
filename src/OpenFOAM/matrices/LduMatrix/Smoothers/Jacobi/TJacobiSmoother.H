#ifndef TJacobiSmoother_H
#define TJacobiSmoother_H

#include "LduMatrix.H"


namespace Foam
{


template<class Type, class DType, class LUType>
class TJacobiSmoother
:
    public LduMatrix<Type, DType, LUType>::smoother
{

    gpuField<DType> rD_;

public:

    TypeName("Jacobi");

    TJacobiSmoother
    (
        const word& fieldName,
        const LduMatrix<Type, DType, LUType>& matrix
    );


    static void smooth
    (
        const word& fieldName,
        gpuField<Type>& psi,
        const LduMatrix<Type, DType, LUType>& matrix,
        const gpuField<DType>& rD,
        const label nSweeps
    );


    virtual void smooth
    (
        gpuField<Type>& psi,
        const label nSweeps
    ) const;
};


} 


#ifdef NoRepository
#   include "TJacobiSmoother.C"
#endif


#endif

