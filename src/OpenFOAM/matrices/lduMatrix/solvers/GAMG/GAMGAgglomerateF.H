#pragma once

namespace Foam
{
namespace GAMG
{


struct negative
{
    __host__ __device__
    bool operator()(const label& x)
    {
        return x<0;
    }
};


struct nonNegative
{
    __host__ __device__
    bool operator()(const label& x)
    {
        return x >= 0;
    }
};


struct faceToDiag
   : public std::unary_function<label,label>
{
    __host__ __device__
    label operator()(const label& n)
    {
        return -1 - n;
    }
};


}
}

