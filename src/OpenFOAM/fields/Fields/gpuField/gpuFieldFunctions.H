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

\*---------------------------------------------------------------------------*/

#define TEMPLATE template<class Type>
#include "gpuFieldFunctionsM.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void component
(
    gpuField<typename gpuField<Type>::cmptType>& res,
    const gpuList<Type>& f,
    const direction d
);


template<class Type>
void T(gpuField<Type>& res, const gpuList<Type>& f);


template<class Type, int r>
void pow
(
    gpuField<typename powProduct<Type, r>::type>& res,
    const gpuList<Type>& vf
);


template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const gpuList<Type>& f,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const tmp<gpuField<Type> >& tf,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);


template<class Type>
void sqr
(
    gpuField<typename outerProduct<Type, Type>::type>& res,
    const gpuList<Type>& vf
);

template<class Type>
tmp<gpuField<typename outerProduct<Type, Type>::type> >
sqr(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename outerProduct<Type, Type>::type> >
sqr(const tmp<gpuField<Type> >& tf);


template<class Type>
void magSqr(gpuField<scalar>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<scalar> > magSqr(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<scalar> > magSqr(const tmp<gpuField<Type> >& tf);


template<class Type>
void mag(gpuField<scalar>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<scalar> > mag(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<scalar> > mag(const tmp<gpuField<Type> >& tf);


template<class Type>
void cmptMax(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMax(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> >
cmptMax(const tmp<gpuField<Type> >& tf);


template<class Type>
void cmptMin(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMin(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> >
cmptMin(const tmp<gpuField<Type> >& tf);


template<class Type>
void cmptAv(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptAv(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptAv(const tmp<gpuField<Type> >& tf);


template<class Type>
void cmptMag(gpuField<Type>& res, const gpuList<Type>& f);

template<class Type>
tmp<gpuField<Type> > cmptMag(const gpuList<Type>& f);

template<class Type>
tmp<gpuField<Type> > cmptMag(const tmp<gpuField<Type> >& tf);

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
                                                                              \
template<class Type>                                                          \
ReturnType Func(const tmp<gpuField<Type> >& tf1);

template<class Type>
Type max(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, sum)

template<class Type>
Type maxMagSqr(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, minMagSqr)


template<class Type>
scalar sumProd(const gpuList<Type>& f1, const gpuList<Type>& f2);

template<class Type>
Type sumCmptProd(const gpuList<Type>& f1, const gpuList<Type>& f2);

template<class Type>
scalar sumSqr(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

template<class Type>
Type sumCmptMag(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const gpuList<Type>& f);

TMP_UNARY_FUNCTION(Type, average)

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                        \
                                                                                \
template<class Type>                                                            \
ReturnType gFunc(const gpuList<Type>& f, const int comm = UPstream::worldComm); \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
scalar gSumProd
(
    const gpuList<Type>& f1,
    const gpuList<Type>& f2,
    const int comm = UPstream::worldComm
);

template<class Type>
Type gSumCmptProd
(
    const gpuList<Type>& f1,
    const gpuList<Type>& f2,
    const int comm = UPstream::worldComm
);

template<class Type>
Type gAverage
(
    const gpuList<Type>& f,
    const int comm = UPstream::worldComm
);

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                     \
                                                                                  \
template<class Type1, class Type2>                                                \
void OpFunc                                                                       \
(                                                                                 \
    gpuField<typename product<Type1, Type2>::type>& res,                          \
    const gpuList<Type1>& f1,                                                     \
    const gpuList<Type2>& f2                                                      \
);                                                                                \
                                                                                  \
template<class Type1, class Type2>                                                \
tmp<gpuField<typename product<Type1, Type2>::type> >                              \
operator Op(const gpuList<Type1>& f1, const gpuList<Type2>& f2);                  \
                                                                                  \
template<class Type1, class Type2>                                                \
tmp<gpuField<typename product<Type1, Type2>::type> >                              \
operator Op(const gpuList<Type1>& f1, const tmp<gpuField<Type2> >& tf2);          \
                                                                                  \
template<class Type1, class Type2>                                                \
tmp<gpuField<typename product<Type1, Type2>::type> >                              \
operator Op(const tmp<gpuField<Type1> >& tf1, const gpuList<Type2>& f2);          \
                                                                                  \
template<class Type1, class Type2>                                                \
tmp<gpuField<typename product<Type1, Type2>::type> >                              \
operator Op(const tmp<gpuField<Type1> >& tf1, const tmp<gpuField<Type2> >& tf2);  \
                                                                                  \
template<class Type, class Form, class Cmpt, int nCmpt>                           \
void OpFunc                                                                       \
(                                                                                 \
    gpuField<typename product<Type, Form>::type>& res,                            \
    const gpuList<Type>& f1,                                                      \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                        \
);                                                                                \
                                                                                  \
template<class Type, class Form, class Cmpt, int nCmpt>                           \
tmp<gpuField<typename product<Type, Form>::type> >                                \
operator Op(const gpuList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs);     \
                                                                                  \
template<class Type, class Form, class Cmpt, int nCmpt>                           \
tmp<gpuField<typename product<Type, Form>::type> >                                \
operator Op(const tmp<gpuField<Type> >&tf1,const VectorSpace<Form,Cmpt,nCmpt>&vs);\
                                                                                  \
template<class Form, class Cmpt, int nCmpt, class Type>                           \
void OpFunc                                                                       \
(                                                                                 \
    gpuField<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                       \
    const gpuList<Type>& f1                                                       \
);                                                                                \
                                                                                  \
template<class Form, class Cmpt, int nCmpt, class Type>                           \
tmp<gpuField<typename product<Form, Type>::type> >                                \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const gpuList<Type>& f1);     \
                                                                                  \
template<class Form, class Cmpt, int nCmpt, class Type>                           \
tmp<gpuField<typename product<Form, Type>::type> >                                \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&vs,const tmp<gpuField<Type> >&tf1);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"
#include "scalarField.H"

// ************************************************************************* //
