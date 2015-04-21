#pragma once

#include "gpuField.H"
#include "BasicCache.H"

namespace Foam
{

class fvMatrixCache
{
    static PtrList<scalargpuField> diagCache;
    static PtrList<scalargpuField> sourceCache;
    static PtrList<scalargpuField> psiCache;

    public:

    static const scalargpuField& diag(label level, label size)
    {
        return cache::retrieve(diagCache,level,size);
    }

    static const scalargpuField& source(label level, label size)
    {
        return cache::retrieve(sourceCache,level,size);
    }

    static const scalargpuField& psi(label level, label size)
    {
        return cache::retrieve(psiCache,level,size);
    }
};

}
