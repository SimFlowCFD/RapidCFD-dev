#pragma once

#include "scalarField.H"

namespace Foam
{


class lduMatrixSolutionCache
{
    static scalargpuField first_;
    static scalargpuField second_;

    static void ensureSize(label size,scalargpuField& field)
    {
        if(field.size() < size)
            field.setSize(size);
    }

public:

    static label favourSpeed;

    static const scalargpuField& first(label size)
    {
        ensureSize(size,first_);
        return first_;
    }

    static const scalargpuField& second(label size)
    {
        ensureSize(size,second_);
        return second_;
    }
};


}


