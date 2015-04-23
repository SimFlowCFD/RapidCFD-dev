#include "JacobiCache.H"

namespace Foam
{
    PtrList<scalargpuField> JacobiCache::psiCache(1);
    PtrList<scalargpuField> JacobiCache::sourceCache(1);
}
