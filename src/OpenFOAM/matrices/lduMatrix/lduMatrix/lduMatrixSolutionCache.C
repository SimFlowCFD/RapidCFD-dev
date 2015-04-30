#include "lduMatrixSolutionCache.H"
#include "debug.H"

namespace Foam
{
    label lduMatrixSolutionCache::favourSpeed
    (
        debug::optimisationSwitch("favourSpeedOverMemory")
    );

    scalargpuField lduMatrixSolutionCache::first_(0);
    scalargpuField lduMatrixSolutionCache::second_(0);
}
