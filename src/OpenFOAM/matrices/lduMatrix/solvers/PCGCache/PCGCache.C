#include "PCGCache.H"

namespace Foam
{
    PtrList<scalargpuField> PCGCache::pACache(1);
    PtrList<scalargpuField> PCGCache::wACache(1);
    PtrList<scalargpuField> PCGCache::rACache(1);

    PtrList<scalargpuField> PCGCache::pTCache(1);
    PtrList<scalargpuField> PCGCache::wTCache(1);
    PtrList<scalargpuField> PCGCache::rTCache(1);
}
