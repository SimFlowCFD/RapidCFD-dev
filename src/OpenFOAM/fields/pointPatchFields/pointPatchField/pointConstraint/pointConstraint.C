#include "pointConstraint.H"
#include "gpuList.C"


namespace Foam
{

template class gpuList<pointConstraint>;

const char* const pTraits<pointConstraint>::typeName = "pointConstraint";
const pointConstraint pTraits<pointConstraint>::zero = pointConstraint();

}

