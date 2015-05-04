#include "fvMatrixCache.H"

namespace Foam
{
    PtrList<scalargpuField> fvMatrixCache::firstCache(1);
    PtrList<scalargpuField> fvMatrixCache::secondCache(1);
    PtrList<scalargpuField> fvMatrixCache::thirdCache(1);
}
