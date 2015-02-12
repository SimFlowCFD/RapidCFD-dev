#include "TJacobiSmoother.H"

namespace Foam
{

struct LduMatrixJacobiFunctor
{
    template<class Type,class Tuple>
    __HOST____DEVICE__
    Type operator()(const Type& psi, const Tuple& t)
    {
        return psi - dot(thrust::get<0>(t),thrust::get<1>(t) - thrust::get<2>(t));
    }
};

}




template<class Type, class DType, class LUType>
Foam::TJacobiSmoother<Type, DType, LUType>::TJacobiSmoother
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix
)
:
    LduMatrix<Type, DType, LUType>::smoother
    (
        fieldName,
        matrix
    ),
    rD_(matrix.diag().size())
{
    const gpuField<DType>& diag = matrix.diag();

    thrust::transform
    (
        diag.begin(),
        diag.end(),
        rD_.begin(),
        invUnaryFunctionFunctor<DType,DType>()
    );
}



template<class Type, class DType, class LUType>
void Foam::TJacobiSmoother<Type, DType, LUType>::smooth
(
    const word& fieldName_,
    gpuField<Type>& psi,
    const LduMatrix<Type, DType, LUType>& matrix_,
    const gpuField<DType>& rD_,
    const label nSweeps
)
{
    //parallel boundaries are not treated correctly
    notImplemented("TJacobiSmoother::smooth");

    // Temporary storage for the product
    gpuField<Type> Apsi(rD_.size());

    const gpuField<Type>& source = matrix_.source();

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        matrix_.Amul
        (
            Apsi,
            psi
        );

        thrust::transform
        (
            psi.begin(),
            psi.end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                rD_.begin(),
                Apsi.begin(),
                source.begin()
            )),
            psi.begin(),
            LduMatrixJacobiFunctor()
        ); 
    }
}


template<class Type, class DType, class LUType>
void Foam::TJacobiSmoother<Type, DType, LUType>::smooth
(
    gpuField<Type>& psi,
    const label nSweeps
) const
{
    smooth(this->fieldName_, psi, this->matrix_, rD_, nSweeps);
}

