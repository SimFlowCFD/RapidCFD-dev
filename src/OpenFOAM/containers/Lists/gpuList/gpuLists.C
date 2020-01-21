#include "gpuList.C"

namespace Foam
{
    template class gpuList<bool>;
    template class gpuList<char>;
    template class gpuList<label>;
    template class gpuList<float>;
    template class gpuList<double>;
}
