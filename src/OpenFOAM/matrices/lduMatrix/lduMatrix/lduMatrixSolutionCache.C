#include "lduMatrixSolutionCache.H"
#include "debug.H"

namespace Foam
{
    label lduMatrixSolutionCache::favourSpeed
    (
        debug::optimisationSwitch("favourSpeedOverMemory")
    );

    scalargpuField lduMatrixSolutionCache::first_;
    scalargpuField lduMatrixSolutionCache::second_;
}
