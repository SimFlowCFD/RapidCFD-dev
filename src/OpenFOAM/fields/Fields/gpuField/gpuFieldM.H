/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    High performance macro functions for gpuField\<Type\> algebra.  These expand
    using either array element access (for vector machines) or pointer
    dereferencing for scalar machines as appropriate.

\*---------------------------------------------------------------------------*/

#ifndef gpuFieldM_H
#define gpuFieldM_H

#include "error.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FULLDEBUG

template<class Type1, class Type2>
void checkgpuFields
(
    const gpuList<Type1>& f1,
    const gpuList<Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorIn
        (
            "checkgpuFields(const gpuList<Type1>&, "
            "const gpuList<Type2>&, const char*)"
        )   << "    incompatible fields"
            << " gpuField<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and gpuField<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << endl << " for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3>
void checkgpuFields
(
    const gpuList<Type1>& f1,
    const gpuList<Type2>& f2,
    const gpuList<Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorIn
        (
            "checkgpuFields(const gpuList<Type1>&, "
            "const gpuList<Type2>&, const gpuList<Type3>&, "
            "const char*)"
        )   << "    incompatible fields"
            << " gpuField<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", gpuField<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and gpuField<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << endl << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<class Type1, class Type2>
void checkgpuFields
(
    const gpuList<Type1>&,
    const gpuList<Type2>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3>
void checkgpuFields
(
    const gpuList<Type1>&,
    const gpuList<Type2>&,
    const gpuList<Type3>&,
    const char*
)
{}

#endif

#define GENERATE_UNARY_OPERATOR_FUNCTORS(Op,opFunc)                                                 \
                                                                                                    \
template<class Type,class RType>                                                                    \
struct opFunc##UnaryOperatorFunctor : public thrust::unary_function<Type, RType>{                   \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type& t){                                                                \
        return Op t;                                                                                \
    }                                                                                               \
};                                                                                                  \

#define GENERATE_OPERATOR_FUNCTORS(Op,opFunc)                                                       \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##OperatorFunctor : public thrust::binary_function<Type1,Type2,RType>{                 \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type1& t1, const Type2& t2){                                             \
        return t1 Op t2;                                                                            \
    }                                                                                               \
};                                                                                                  \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##OperatorFSFunctor : public thrust::unary_function<Type1, RType>{                     \
    const Type2 t2;                                                                                 \
    opFunc##OperatorFSFunctor(Type2 _t2):t2(_t2){}                                                  \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type1& t1){                                                              \
        return t1 Op t2;                                                                            \
    }                                                                                               \
};                                                                                                  \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##OperatorSFFunctor : public thrust::unary_function<Type2, RType>{                     \
    const Type1 t1;                                                                                 \
    opFunc##OperatorSFFunctor(Type1 _t1):t1(_t1){}                                                  \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type2& t2){                                                              \
        return t1 Op t2;                                                                            \
    }                                                                                               \
};                                                                                                  \
                                                                                                    \

#define GENERATE_UNARY_FUNCTION_FUNCTORS(opFunc)                                                    \
                                                                                                    \
template<class Type,class RType>                                                                    \
struct opFunc##UnaryFunctionFunctor : public thrust::unary_function<Type, RType>{                   \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type& t) const {                                                         \
        return opFunc(t);                                                                           \
    }                                                                                               \
};                                                                                                  \

#define GENERATE_BINARY_FUNCTION_FUNCTORS(opFunc)                                                   \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##BinaryFunctionFunctor : public thrust::binary_function<Type1,Type2,RType>{           \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type1& t1, const Type2& t2) const {                                      \
        return opFunc(t1, t2);                                                                      \
    }                                                                                               \
};                                                                                                  \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##BinaryFunctionFSFunctor : public thrust::unary_function<Type1, RType>{               \
    const Type2 t2;                                                                                 \
    opFunc##BinaryFunctionFSFunctor(Type2 _t2):t2(_t2){}                                            \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type1& t1) const {                                                       \
        return opFunc(t1, t2);                                                                      \
    }                                                                                               \
};                                                                                                  \
                                                                                                    \
template<class Type1,class Type2,class RType>                                                       \
struct opFunc##BinaryFunctionSFFunctor : public thrust::unary_function<Type2, RType>{               \
    const Type1 t1;                                                                                 \
    opFunc##BinaryFunctionSFFunctor(Type1 _t1):t1(_t1){}                                            \
    __HOST____DEVICE__                                                                              \
    RType operator()(const Type2& t2) const {                                                       \
        return opFunc(t1, t2);                                                                      \
    }                                                                                               \
};                                                                                                  \
 
GENERATE_UNARY_OPERATOR_FUNCTORS(-, negate)
GENERATE_UNARY_OPERATOR_FUNCTORS(*, hdual)

GENERATE_OPERATOR_FUNCTORS(+,      add)
GENERATE_OPERATOR_FUNCTORS(-, subtract)
GENERATE_OPERATOR_FUNCTORS(*, multiply)
GENERATE_OPERATOR_FUNCTORS(/,   divide)
GENERATE_OPERATOR_FUNCTORS(*,    outer)
GENERATE_OPERATOR_FUNCTORS(^,    cross)
GENERATE_OPERATOR_FUNCTORS(&,      dot)
GENERATE_OPERATOR_FUNCTORS(&&,  dotdot)

GENERATE_BINARY_FUNCTION_FUNCTORS(dot)
GENERATE_BINARY_FUNCTION_FUNCTORS(min)
GENERATE_BINARY_FUNCTION_FUNCTORS(minMod)
GENERATE_BINARY_FUNCTION_FUNCTORS(max)
GENERATE_BINARY_FUNCTION_FUNCTORS(pow)
GENERATE_BINARY_FUNCTION_FUNCTORS(atan2)
GENERATE_BINARY_FUNCTION_FUNCTORS(cmptMultiply)
GENERATE_BINARY_FUNCTION_FUNCTORS(cmptDivide)
GENERATE_BINARY_FUNCTION_FUNCTORS(transform)


GENERATE_UNARY_FUNCTION_FUNCTORS(diag)
GENERATE_UNARY_FUNCTION_FUNCTORS(pow3)
GENERATE_UNARY_FUNCTION_FUNCTORS(pow4)
GENERATE_UNARY_FUNCTION_FUNCTORS(pow5)
GENERATE_UNARY_FUNCTION_FUNCTORS(pow6)
GENERATE_UNARY_FUNCTION_FUNCTORS(pow025)
GENERATE_UNARY_FUNCTION_FUNCTORS(sqrt)
GENERATE_UNARY_FUNCTION_FUNCTORS(cbrt)
GENERATE_UNARY_FUNCTION_FUNCTORS(sign)
GENERATE_UNARY_FUNCTION_FUNCTORS(pos)
GENERATE_UNARY_FUNCTION_FUNCTORS(neg)
GENERATE_UNARY_FUNCTION_FUNCTORS(exp)
GENERATE_UNARY_FUNCTION_FUNCTORS(log)
GENERATE_UNARY_FUNCTION_FUNCTORS(log10)
GENERATE_UNARY_FUNCTION_FUNCTORS(sin)
GENERATE_UNARY_FUNCTION_FUNCTORS(cos)
GENERATE_UNARY_FUNCTION_FUNCTORS(tan)
GENERATE_UNARY_FUNCTION_FUNCTORS(asin)
GENERATE_UNARY_FUNCTION_FUNCTORS(acos)
GENERATE_UNARY_FUNCTION_FUNCTORS(atan)
GENERATE_UNARY_FUNCTION_FUNCTORS(sinh)
GENERATE_UNARY_FUNCTION_FUNCTORS(cosh)
GENERATE_UNARY_FUNCTION_FUNCTORS(tanh)
GENERATE_UNARY_FUNCTION_FUNCTORS(asinh)
GENERATE_UNARY_FUNCTION_FUNCTORS(acosh)
GENERATE_UNARY_FUNCTION_FUNCTORS(atanh)
GENERATE_UNARY_FUNCTION_FUNCTORS(erf)
GENERATE_UNARY_FUNCTION_FUNCTORS(erfc)
GENERATE_UNARY_FUNCTION_FUNCTORS(lgamma)
GENERATE_UNARY_FUNCTION_FUNCTORS(j0)
GENERATE_UNARY_FUNCTION_FUNCTORS(j1)
GENERATE_UNARY_FUNCTION_FUNCTORS(y0)
GENERATE_UNARY_FUNCTION_FUNCTORS(y1)

GENERATE_UNARY_FUNCTION_FUNCTORS(degToRad)
GENERATE_UNARY_FUNCTION_FUNCTORS(radToDeg)
GENERATE_UNARY_FUNCTION_FUNCTORS(atmToPa)
GENERATE_UNARY_FUNCTION_FUNCTORS(paToAtm)

GENERATE_UNARY_FUNCTION_FUNCTORS(mag)
GENERATE_UNARY_FUNCTION_FUNCTORS(magSqr)
GENERATE_UNARY_FUNCTION_FUNCTORS(sqr)
GENERATE_UNARY_FUNCTION_FUNCTORS(cmptMag)
GENERATE_UNARY_FUNCTION_FUNCTORS(cmptMax)
GENERATE_UNARY_FUNCTION_FUNCTORS(cmptMin)
GENERATE_UNARY_FUNCTION_FUNCTORS(cmptAv)
GENERATE_UNARY_FUNCTION_FUNCTORS(eigenVectors)
GENERATE_UNARY_FUNCTION_FUNCTORS(eigenValues)

GENERATE_UNARY_FUNCTION_FUNCTORS(tr)
GENERATE_UNARY_FUNCTION_FUNCTORS(sph)
GENERATE_UNARY_FUNCTION_FUNCTORS(symm)
GENERATE_UNARY_FUNCTION_FUNCTORS(twoSymm)
GENERATE_UNARY_FUNCTION_FUNCTORS(skew)
GENERATE_UNARY_FUNCTION_FUNCTORS(dev)
GENERATE_UNARY_FUNCTION_FUNCTORS(dev2)
GENERATE_UNARY_FUNCTION_FUNCTORS(det)
GENERATE_UNARY_FUNCTION_FUNCTORS(cof)
GENERATE_UNARY_FUNCTION_FUNCTORS(inv)


#undef GENERATE_UNARY_OPERATOR_FUNCTORS
#undef GENERATE_OPERATOR_FUNCTORS
#undef GENERATE_BINARY_FUNCTION_FUNCTORS
#undef GENERATE_UNARY_FUNCTION_FUNCTORS

template<class Type,class RType>
struct assignFunctor
{
    __host__ __device__
    RType operator()(const Type& t) const 
    {
        return t;
    }
};

template<class Type>
struct transposeFunctor
{
    __host__ __device__
    Type operator()(const Type& t) const 
    {
        return t.T();
    }
};

struct stabiliseFunctor
{
    const scalar sm;
    stabiliseFunctor(scalar _s): sm(_s){}
    __host__ __device__
    scalar operator()(const scalar& sf)
    {
        if (sf >= 0)
        {
            return sf + sm;
        }
        else
        {
            return sf - sm;
        }
    }
};

template<class Type1,class Type2>
struct assignValueFunctor
{
    const Type2 t2;
    assignValueFunctor(Type2 _t2): t2(_t2) {}
    __host__ __device__
    Type1 operator()() const 
    {
        return t2;
    }
};

template<class Type1,class Type2>
struct componentFunctor
{
    const direction d;
    componentFunctor(direction _d): d(_d) {}
    __host__ __device__
    Type1 operator()(const Type2& tt)
    {
        return tt.component(d);
    }
};

template<class Type,class cmptType>
struct replaceComponentFunctor
{
    const direction d;
    replaceComponentFunctor(direction _d): d(_d) {}
    __host__ __device__
    Type operator()(const Type& tt,const cmptType& ct)
    {
        Type tr(tt);
        tr.replace(d,ct);
        return tr;
    }
};

template<class Type,class cmptType>
struct replaceComponentWithSourceFunctor
{
    const direction d;
    const cmptType ct;
    replaceComponentWithSourceFunctor
    (
        direction _d,
        cmptType _ct
    ): 
        d(_d), 
        ct(_ct) 
    {}

    __host__ __device__
    Type operator()(const Type& tt)
    {
        Type tr(tt);
        tr.replace(d,ct);
        return tr;
    }
};

template<class Type>
struct addWithSourceFunctor
{
    const Type t;
    addWithSourceFunctor(Type _t): t(_t) {}
    __host__ __device__
    Type operator()(const Type& t1)
    {
        return t + t1;
    }
};

template<class Type1,class Type2>
struct powerFunctor
{
    const Type2 tt;
    powerFunctor(Type2 _tt):tt(_tt){}
    __host__ __device__
    Type2 operator()(const Type1& t)
    {
        return pow(t,tt);
    }
};

template<class Type1,class Type2>
struct outerProductFunctor
{
    __host__ __device__
    Type2 operator()(const Type1& t)
    {
        return sqr(t);
    }
};

template<class Type>
struct productFunctor
{
    __host__ __device__
    scalar operator()(const Type& t1,const Type& t2)
    {
        return t1 && t2;
    }
};

template<class Type>
struct maximumMagnitudeSquaredFunctor
{
    __host__ __device__
    Type operator()(const Type& t1,const Type& t2)
    {
        return maxMagSqrOp<Type>()(t1,t2);
    }
};

template<class Type>
struct minimumMagnitudeSquaredFunctor
{
    __host__ __device__
    Type operator()(const Type& t1,const Type& t2)
    {
        return minMagSqrOp<Type>()(t1,t2);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
