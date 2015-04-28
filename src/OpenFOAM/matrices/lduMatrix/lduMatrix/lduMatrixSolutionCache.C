#include "lduMatrixSolutionCache.H"

namespace Foam
{
    scalargpuField lduMatrixSolutionCache::first_(0);
    scalargpuField lduMatrixSolutionCache::second_(0);
}
