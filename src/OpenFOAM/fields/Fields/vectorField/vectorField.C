#include "vectorField.H"
#include "gpuList.C"

#define TEMPLATE
#include "gpuFieldFunctionsM.C"

namespace Foam
{

template class gpuList<vector>;
template class gpuField<vector>;

template Ostream& operator<<<vector>(Ostream&, const gpuList<vector>&);
template Istream& operator>><vector>(Istream&, gpuList<vector>&);

template<>
__host__ __device__
vector transposeFunctor<vector>::operator()(const vector& v) const
{
    return v;
}

BINARY_SYM_OPERATOR(vector, scalar, vector, *, outer)
BINARY_SYM_FUNCTION(vector, scalar, vector, multiply)
BINARY_OPERATOR(vector, vector, scalar, /, divide)
BINARY_TYPE_OPERATOR_FS(vector, vector, scalar, /, divide)

BINARY_FULL_OPERATOR(vector, vector, vector, +, add)
BINARY_FULL_OPERATOR(vector, vector, vector, -, subtract)

BINARY_FULL_OPERATOR(scalar, vector, vector, &, dot)
BINARY_FULL_OPERATOR(vector, vector, vector, ^, cross)

}


#include "undefgpuFieldFunctionsM.H"

#include "gpuFieldCommonFunctions.C"
// force instantiation
#define TEMPLATE template
#define FTYPE vector
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"
