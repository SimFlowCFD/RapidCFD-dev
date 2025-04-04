#include "PCGCache.H"

namespace Foam
{
    PtrList<scalargpuField> PCGCache::pACache(1);
    PtrList<scalargpuField> PCGCache::wACache(1);
    PtrList<scalargpuField> PCGCache::rACache(1);

    PtrList<scalargpuField> PCGCache::pTCache(1);
    PtrList<scalargpuField> PCGCache::wTCache(1);
    PtrList<scalargpuField> PCGCache::rTCache(1);
    /*
    PtrList<scalargpuField> PCGCache::yACache(1);
    PtrList<scalargpuField> PCGCache::AyACache(1);
    PtrList<scalargpuField> PCGCache::sACache(1);
    PtrList<scalargpuField> PCGCache::zACache(1);
    */
    PtrList<scalargpuField> PCGCache::tACache(1);
    PtrList<scalargpuField> PCGCache::result1Cache(1);
}
