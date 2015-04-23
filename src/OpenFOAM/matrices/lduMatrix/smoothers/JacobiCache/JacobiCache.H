#pragma once

#include "BasicCache.H"

namespace Foam
{
    class JacobiCache
    {
        static PtrList<scalargpuField> psiCache;
        static PtrList<scalargpuField> sourceCache;

        public:

        static const scalargpuField& psi(label level, label size)
        {
            return cache::retrieve(psiCache,level,size);
        }

        static const scalargpuField& source(label level, label size)
        {
            return cache::retrieve(sourceCache,level,size);
        }
    };
}


