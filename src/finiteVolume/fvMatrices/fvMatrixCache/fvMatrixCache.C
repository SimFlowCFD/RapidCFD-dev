#include "fvMatrixCache.H"

namespace Foam
{
    PtrList<scalargpuField> fvMatrixCache::diagCache(1);
    PtrList<scalargpuField> fvMatrixCache::sourceCache(1);
    PtrList<scalargpuField> fvMatrixCache::psiCache(1);
}
