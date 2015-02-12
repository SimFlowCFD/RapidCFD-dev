/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.

    The information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define EqOp(opName, op)                                                    \
                                                                            \
template<class T1, class T2>                                                \
class opName##Op2                                                           \
{                                                                           \
public:                                                                     \
                                                                            \
    __HOST____DEVICE__                                                     \
    void operator()(T1& x, const T2& y) const                               \
    {                                                                       \
        op;                                                                 \
    }                                                                       \
};                                                                          \
                                                                            \
template<class T>                                                           \
class opName##Op                                                            \
{                                                                           \
public:                                                                     \
                                                                            \
    __HOST____DEVICE__                                                     \
    void operator()(T& x, const T& y) const                                 \
    {                                                                       \
        op;                                                                 \
    }                                                                       \
};

EqOp(eq, x = y)
EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(eqMag, x = mag(y))
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(maxEq, x = max(x, y))
EqOp(minEq, x = min(x, y))
EqOp(minMagSqrEq, x = (magSqr(x)<=magSqr(y) ? x : y))
EqOp(maxMagSqrEq, x = (magSqr(x)>=magSqr(y) ? x : y))
EqOp(andEq, x = (x && y))
EqOp(orEq, x  = (x || y))

EqOp(eqMinus, x = -y)

EqOp(nopEq, (void)x)

#undef EqOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if __GNUC__
#define WARNRETURN __attribute__((warn_unused_result))
#else
#define WARNRETURN
#endif

#define Op(opName, op)                                                        \
                                                                              \
    template<class T, class T1, class T2>                                     \
    class opName##Op3                                                         \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T operator()(const T1& x, const T2& y) const WARNRETURN               \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template<class T1, class T2>                                              \
    class opName##Op2                                                         \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T1 operator()(const T1& x, const T2& y) const WARNRETURN              \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template<class T>                                                         \
    class opName##Op                                                          \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T operator()(const T& x, const T& y) const WARNRETURN                 \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template<class T, class T1, class T2>                                     \
    class opName##SwapOp3                                                     \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T operator()(const T1& y, const T2& x) const WARNRETURN               \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template<class T1, class T2>                                              \
    class opName##SwapOp2                                                     \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T1 operator()(const T1& y, const T2& x) const WARNRETURN              \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template<class T>                                                         \
    class opName##SwapOp                                                      \
    {                                                                         \
    public:                                                                   \
                                                                              \
        __HOST____DEVICE__                                                    \
        T operator()(const T& y, const T& x) const WARNRETURN                 \
        {                                                                     \
            return op;                                                        \
        }                                                                     \
    };

#define weightedOp(opName, op)                                                \
                                                                              \
    template<class Type, class CombineOp>                                     \
    class opName##WeightedOp                                                  \
    {                                                                         \
        const CombineOp& cop_;                                                \
                                                                              \
        public:                                                               \
                                                                              \
            opName##WeightedOp(const CombineOp& cop)                          \
            :                                                                 \
                cop_(cop)                                                     \
            {}                                                                \
                                                                              \
            __HOST____DEVICE__                                                \
            void operator()                                                   \
            (                                                                 \
                Type& x,                                                      \
                const label index,                                            \
                const Type& y,                                                \
                const scalar weight                                           \
            ) const                                                           \
            {                                                                 \
                cop_(x, op);                                                  \
            }                                                                 \
    };                                                                        \


Op(sum, x + y)
Op(right, y)
Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(max, max(x, y))
Op(min, min(x, y))
Op(minMagSqr, (magSqr(x)<=magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x)>=magSqr(y) ? x : y))
Op(minMod, minMod(x, y))
Op(and, x && y)
Op(or, x || y)
Op(eqEq, x == y)
Op(less, x < y)
Op(lessEq, x <= y)
Op(greater, x > y)
Op(greaterEq, x >= y)

weightedOp(multiply, (weight*y))

template<class T>                                              
class unityOp                                                         
{                                                                         
public:                                                                   
                                                                              
    __HOST____DEVICE__                                                    
    T operator()(const T& x) const            
    {                                                                     
        return x;                                                        
    }                                                                     
 };  

template<class T>                                              
class magSqrOp                                                         
{                                                                         
public:                                                                       
    __HOST____DEVICE__                                                    
    T operator()(const T& x) const            
    {                                                                     
        return magSqr(x);                                                        
    }                                                                     
};   

template<class T>                                              
class magOp                                                         
{                                                                         
public:                                                                       
    __HOST____DEVICE__                                                    
    T operator()(const T& x) const            
    {                                                                     
        return mag(x);                                                        
    }                                                                     
};                                                                      

#undef Op
#undef weightedOp
#undef WARNRETURN

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
