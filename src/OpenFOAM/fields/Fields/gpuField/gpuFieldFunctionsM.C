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

\*---------------------------------------------------------------------------*/

#include "gpuFieldM.H"
#include "gpuFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(gpuField<ReturnType>& res, const gpuList<Type>& f)                   \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f.begin(),                                                             \
        f.end(),                                                               \
        res.begin(),                                                           \
        Func##UnaryFunctionFunctor<Type,ReturnType>()                          \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func(const gpuList<Type>& f)                        \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f.size()));       \
    Func(tRes(), f);                                                           \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func(const tmp<gpuField<Type> >& tf)                \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type>::New(tf);  \
    Func(tRes(), tf());                                                        \
    reusegpuTmp<ReturnType, Type>::clear(tf);                                  \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(gpuField<ReturnType>& res, const gpuList<Type>& f)                 \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f.begin(),                                                             \
        f.end(),                                                               \
        res.begin(),                                                           \
        OpFunc##UnaryOperatorFunctor<Type,ReturnType>()                        \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op(const gpuList<Type>& f)                 \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f.size()));       \
    OpFunc(tRes(), f);                                                         \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op(const tmp<gpuField<Type> >& tf)         \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type>::New(tf);  \
    OpFunc(tRes(), tf());                                                      \
    reusegpuTmp<ReturnType, Type>::clear(tf);                                  \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const gpuList<Type1>& f1,                                                  \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f1.begin(),                                                            \
        f1.end(),                                                              \
        f2.begin(),                                                            \
        res.begin(),                                                           \
        Func##BinaryFunctionFunctor<Type1,Type2,ReturnType>()                  \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f1.size()));      \
    Func(tRes(), f1, f2);                                                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type2>::New(tf2);\
    Func(tRes(), f1, tf2());                                                   \
    reusegpuTmp<ReturnType, Type2>::clear(tf2);                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type1>::New(tf1);\
    Func(tRes(), tf1(), f2);                                                   \
    reusegpuTmp<ReturnType, Type1>::clear(tf1);                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes =                                          \
        reusegpuTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);        \
    Func(tRes(), tf1(), tf2());                                                \
    reusegpuTmpTmp<ReturnType, Type1, Type1, Type2>::clear(tf1, tf2);          \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const Type1& s1,                                                           \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f2.begin(),                                                            \
        f2.end(),                                                              \
        res.begin(),                                                           \
        Func##BinaryFunctionSFFunctor<Type1,Type2,ReturnType>(s1)              \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const Type1& s1,                                                           \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f2.size()));      \
    Func(tRes(), s1, f2);                                                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type2>::New(tf2);\
    Func(tRes(), s1, tf2());                                                   \
    reusegpuTmp<ReturnType, Type2>::clear(tf2);                                \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const gpuList<Type1>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f1.begin(),                                                            \
        f1.end(),                                                              \
        res.begin(),                                                           \
        Func##BinaryFunctionFSFunctor<Type1,Type2,ReturnType>(s2)              \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f1.size()));      \
    Func(tRes(), f1, s2);                                                      \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > Func                                                \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type1>::New(tf1);\
    Func(tRes(), tf1(), s2);                                                   \
    reusegpuTmp<ReturnType, Type1>::clear(tf1);                                \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const gpuList<Type1>& f1,                                                  \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f1.begin(),                                                            \
        f1.end(),                                                              \
        f2.begin(),                                                            \
        res.begin(),                                                           \
        OpFunc##OperatorFunctor<Type1,Type2,ReturnType>()                      \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f1.size()));      \
    OpFunc(tRes(), f1, f2);                                                    \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type2>::New(tf2);\
    OpFunc(tRes(), f1, tf2());                                                 \
    reusegpuTmp<ReturnType, Type2>::clear(tf2);                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type1>::New(tf1);\
    OpFunc(tRes(), tf1(), f2);                                                 \
    reusegpuTmp<ReturnType, Type1>::clear(tf1);                                \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes =                                          \
        reusegpuTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);        \
    OpFunc(tRes(), tf1(), tf2());                                              \
    reusegpuTmpTmp<ReturnType, Type1, Type1, Type2>::clear(tf1, tf2);          \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const Type1& s1,                                                           \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f2.begin(),                                                            \
        f2.end(),                                                              \
        res.begin(),                                                           \
        OpFunc##OperatorSFFunctor<Type1,Type2,ReturnType>(s1)                  \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const Type1& s1,                                                           \
    const gpuList<Type2>& f2                                                   \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f2.size()));      \
    OpFunc(tRes(), s1, f2);                                                    \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<gpuField<Type2> >& tf2                                           \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type2>::New(tf2);\
    OpFunc(tRes(), s1, tf2());                                                 \
    reusegpuTmp<ReturnType, Type2>::clear(tf2);                                \
    return tRes;                                                               \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    gpuField<ReturnType>& res,                                                 \
    const gpuList<Type1>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    thrust::transform                                                          \
    (                                                                          \
        f1.begin(),                                                            \
        f1.end(),                                                              \
        res.begin(),                                                           \
        OpFunc##OperatorFSFunctor<Type1,Type2,ReturnType>(s2)                  \
    );                                                                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const gpuList<Type1>& f1,                                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes(new gpuField<ReturnType>(f1.size()));      \
    OpFunc(tRes(), f1, s2);                                                    \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<gpuField<ReturnType> > operator Op                                         \
(                                                                              \
    const tmp<gpuField<Type1> >& tf1,                                          \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    tmp<gpuField<ReturnType> > tRes = reusegpuTmp<ReturnType, Type1>::New(tf1);\
    OpFunc(tRes(), tf1(), s2);                                                 \
    reusegpuTmp<ReturnType, Type1>::clear(tf1);                                \
    return tRes;                                                               \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)


// ************************************************************************* //
