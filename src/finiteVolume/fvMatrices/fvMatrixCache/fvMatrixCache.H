#pragma once

#include "gpuField.H"
#include "BasicCache.H"

namespace Foam
{

class fvMatrixCache
{
    static PtrList<scalargpuField> firstCache;
    static PtrList<scalargpuField> secondCache;
    static PtrList<scalargpuField> thirdCache;

    public:

    static const scalargpuField& first(label level, label size)
    {
        return cache::retrieve(firstCache,level,size);
    }

    static const scalargpuField& second(label level, label size)
    {
        return cache::retrieve(secondCache,level,size);
    }

    static const scalargpuField& third(label level, label size)
    {
        return cache::retrieve(thirdCache,level,size);
    }
};

}
