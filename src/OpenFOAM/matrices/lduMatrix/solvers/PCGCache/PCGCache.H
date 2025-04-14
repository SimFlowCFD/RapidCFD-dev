#pragma once

#include "gpuField.H"
#include "BasicCache.H"

namespace Foam
{

class PCGCache
{
    static PtrList<scalargpuField> pACache;
    static PtrList<scalargpuField> wACache;
    static PtrList<scalargpuField> rACache;

    static PtrList<scalargpuField> pTCache;
    static PtrList<scalargpuField> wTCache;
    static PtrList<scalargpuField> rTCache;

    static PtrList<scalargpuField> tACache;
    static PtrList<scalargpuField> result1Cache;

    public:

    static const scalargpuField& pA(label level, label size)
    {
        return cache::retrieveConst(pACache,level,size);
    }

    static const scalargpuField& wA(label level, label size)
    {
        return cache::retrieveConst(wACache,level,size);
    }

    static const scalargpuField& rA(label level, label size)
    {
        return cache::retrieveConst(rACache,level,size);
    }

    static const scalargpuField& pT(label level, label size)
    {
        return cache::retrieveConst(pTCache,level,size);
    }

    static const scalargpuField& wT(label level, label size)
    {
        return cache::retrieveConst(wTCache,level,size);
    }

    static const scalargpuField& rT(label level, label size)
    {
        return cache::retrieveConst(rTCache,level,size);
    }

    static const scalargpuField& tA(label level, label size)
    {
        return cache::retrieveConst(tACache,level,size);
    }

    static const scalargpuField& result1(label level, label size)
    {
        return cache::retrieveConst(result1Cache,level,size);
    }
};

}
