#include "gpuField.H"
#include "vector.H"

namespace Foam
{

template class gpuField<vector>;

template<>
__host__ __device__
vector transposeFunctor<vector>::operator()(const vector& v) const
{
    return v;
}

}
