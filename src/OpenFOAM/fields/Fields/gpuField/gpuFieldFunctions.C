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

#include "PstreamReduceOps.H"
#include "gpuFieldReuseFunctions.H"

#define TEMPLATE template<class Type>
#include "gpuFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    gpuField<typename gpuField<Type>::cmptType>& res,
    const gpuList<Type>& f,
    const direction d
)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        componentFunctor<cmptType,Type>(d)
    );
}

template<class Type>
void T(gpuField<Type>& res, const gpuList<Type>& f)
{
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        transposeFunctor<Type>()
    );
}

template<class Type, int r>
void pow
(
    gpuField<typename powProduct<Type, r>::type>& res,
    const gpuList<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    thrust::transform
    (
        vf.begin(),
        vf.end(),
        res.begin(),
        powBinaryFunctionFSFunctor<Type,powProductType,powProductType>
        (
            pTraits<powProductType>::zero
        )
    );
}

template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const gpuList<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<gpuField<powProductType> > tRes
    (
        new gpuField<powProductType>(f.size())
    );
    pow<Type, r>(tRes(), f);
    return tRes;
}

template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const tmp<gpuField<Type> >& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<gpuField<powProductType> > tRes = reusegpuTmp<powProductType, Type>::New(tf);
    pow<Type, r>(tRes(), tf());
    reusegpuTmp<powProductType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void sqr
(
    gpuField<typename outerProduct<Type, Type>::type>& res,
    const gpuList<Type>& vf
)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    thrust::transform
    (
        vf.begin(),
        vf.end(),
        res.begin(),
        outerProductFunctor<Type,outerProductType>()
    );
}

template<class Type>
tmp<gpuField<typename outerProduct<Type, Type>::type> >
sqr(const gpuList<Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<gpuField<outerProductType> > tRes
    (
        new gpuField<outerProductType>(f.size())
    );
    sqr(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<typename outerProduct<Type, Type>::type> >
sqr(const tmp<gpuField<Type> >& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<gpuField<outerProductType> > tRes =
        reusegpuTmp<outerProductType, Type>::New(tf);
    sqr(tRes(), tf());
    reusegpuTmp<outerProductType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void magSqr(gpuField<scalar>& res, const gpuList<Type>& f)
{
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        magSqrUnaryFunctionFunctor<Type,scalar>()
    );
}

template<class Type>
tmp<gpuField<scalar> > magSqr(const gpuList<Type>& f)
{
    tmp<gpuField<scalar> > tRes(new gpuField<scalar>(f.size()));
    magSqr(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<scalar> > magSqr(const tmp<gpuField<Type> >& tf)
{
    tmp<gpuField<scalar> > tRes = reusegpuTmp<scalar, Type>::New(tf);
    magSqr(tRes(), tf());
    reusegpuTmp<scalar, Type>::clear(tf);
    return tRes;
}


template<class Type>
void mag(gpuField<scalar>& res, const gpuList<Type>& f)
{
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        magUnaryFunctionFunctor<Type,scalar>()
    );
}

template<class Type>
tmp<gpuField<scalar> > mag(const gpuList<Type>& f)
{
    tmp<gpuField<scalar> > tRes(new gpuField<scalar>(f.size()));
    mag(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<scalar> > mag(const tmp<gpuField<Type> >& tf)
{
    tmp<gpuField<scalar> > tRes = reusegpuTmp<scalar, Type>::New(tf);
    mag(tRes(), tf());
    reusegpuTmp<scalar, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMax(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        cmptMaxUnaryFunctionFunctor<Type,cmptType>()
    );
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMax(const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes(new gpuField<cmptType>(f.size()));
    cmptMax(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMax(const tmp<gpuField<Type> >& tf)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes = reusegpuTmp<cmptType, Type>::New(tf);
    cmptMax(tRes(), tf());
    reusegpuTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMin(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        cmptMaxUnaryFunctionFunctor<Type,cmptType>()
    );
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMin(const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes(new gpuField<cmptType>(f.size()));
    cmptMin(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptMin(const tmp<gpuField<Type> >& tf)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes = reusegpuTmp<cmptType, Type>::New(tf);
    cmptMin(tRes(), tf());
    reusegpuTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptAv(gpuField<typename gpuField<Type>::cmptType>& res, const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        cmptAvUnaryFunctionFunctor<Type,cmptType>()
    );
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptAv(const gpuList<Type>& f)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes(new gpuField<cmptType>(f.size()));
    cmptAv(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<typename gpuField<Type>::cmptType> > cmptAv(const tmp<gpuField<Type> >& tf)
{
    typedef typename gpuField<Type>::cmptType cmptType;
    tmp<gpuField<cmptType> > tRes = reusegpuTmp<cmptType, Type>::New(tf);
    cmptAv(tRes(), tf());
    reusegpuTmp<cmptType, Type>::clear(tf);
    return tRes;
}


template<class Type>
void cmptMag(gpuField<Type>& res, const gpuList<Type>& f)
{
    thrust::transform
    (
        f.begin(),
        f.end(),
        res.begin(),
        cmptMagUnaryFunctionFunctor<Type,Type>()
    );
}

template<class Type>
tmp<gpuField<Type> > cmptMag(const gpuList<Type>& f)
{
    tmp<gpuField<Type> > tRes(new gpuField<Type>(f.size()));
    cmptMag(tRes(), f);
    return tRes;
}

template<class Type>
tmp<gpuField<Type> > cmptMag(const tmp<gpuField<Type> >& tf)
{
    tmp<gpuField<Type> > tRes = reusegpuTmp<Type, Type>::New(tf);
    cmptMag(tRes(), tf());
    reusegpuTmp<Type, Type>::clear(tf);
    return tRes;
}


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
                                                                              \
template<class Type>                                                          \
ReturnType Func(const tmp<gpuField<Type> >& tf1)                              \
{                                                                             \
    ReturnType res = Func(tf1());                                             \
    tf1.clear();                                                              \
    return res;                                                               \
}

template<class Type>
Type max(const gpuList<Type>& f)
{
    if (f.size())
    {
        return *thrust::max_element(f.begin(), f.end());
    }
    else
    {
        return pTraits<Type>::min;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const gpuList<Type>& f)
{
    if (f.size())
    {
        return *thrust::min_element(f.begin(), f.end());
    }
    else
    {
        return pTraits<Type>::max;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const gpuList<Type>& f)
{
    if (f.size())
    {
        return thrust::reduce
               (
                   f.begin(), 
                   f.end(), 
                   pTraits<Type>::zero, 
                   thrust::plus<Type>()
               );
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)

template<class Type>
Type maxMagSqr(const gpuList<Type>& f)
{
    if (f.size())
    {
        gpuList<scalar>& ms(f.size());

	thrust::transform
        (
            f.begin(),
            f.end(),
            ms.begin(),
            magSqrUnaryFunctionFunctor<Type,scalar>()
        );
	typename thrust::device_vector<scalar>::iterator iter = 
                     thrust::max_element(ms.begin(), ms.end());

        unsigned int position = iter - ms.begin();

        return f.get(position);
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

template<class Type>
Type minMagSqr(const gpuList<Type>& f)
{
    if (f.size())
    {
        gpuList<scalar>& ms (f.size());

	thrust::transform
        (
            f.begin(),
            f.end(),
            ms.begin(),
            magSqrUnaryFunctionFunctor<Type,scalar>()
        );

	typename thrust::device_vector<scalar>::iterator  iter = 
                     thrust::min_element(ms.begin(), ms.end());

        unsigned int position = iter - ms.begin();

        return f.get(position);
    }
    else
    {
        return pTraits<Type>::rootMax;
    }
}

template<class Type>
scalar sumProd(const gpuList<Type>& f1, const gpuList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        gpuList<scalar> tmp(f1.size());

        thrust::transform
        (
            f1.begin(),
            f1.end(),f2.begin(),
            tmp.begin(),productFunctor<Type>()
        );

        return thrust::reduce
               (
                   tmp.begin(),
                   tmp.end(),
                   pTraits<scalar>::zero,
                   thrust::plus<scalar>()
               );
    }
    else
    {
        return 0.0;
    }
}


template<class Type>
Type sumCmptProd(const gpuList<Type>& f1, const gpuList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    { 
        gpuList<Type> tmp(f1.size());

        thrust::transform
        (
            f1.begin(),
            f1.end(),
            f2.begin(),
            tmp.begin(),
            cmptMultiplyBinaryFunctionFunctor<Type,Type,Type>()
        );

        return thrust::reduce
               (
                   tmp.begin(),
                   tmp.end(),
                   pTraits<Type>::zero,
                   thrust::plus<Type>()
               );
    }
    else
    {
        return pTraits<Type>::zero;
    }
}


template<class Type>
scalar sumSqr(const gpuList<Type>& f)
{
    if (f.size())
    {
        return thrust::reduce
               (
                   thrust::make_transform_iterator
                   (
                       f.begin(),
                       outerProductFunctor<Type,scalar>()
                   ),
                   thrust::make_transform_iterator
                   (
                       f.end(),
                       outerProductFunctor<Type,scalar>()
                   ),
                   pTraits<scalar>::zero,
                   thrust::plus<scalar>()
                );
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const gpuList<Type>& f)
{
    if (f.size())
    {
        return thrust::reduce
               (
                   thrust::make_transform_iterator
                   (
                       f.begin(),
                       magUnaryFunctionFunctor<Type,scalar>()
                   ),
                   thrust::make_transform_iterator
                   (
                       f.end(),
                       magUnaryFunctionFunctor<Type,scalar>()
                   ),
                   pTraits<scalar>::zero,
                   thrust::plus<scalar>()
               );
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumMag)


template<class Type>
Type sumCmptMag(const gpuList<Type>& f)
{
    if (f.size())
    {
        return thrust::reduce
               (
                   thrust::make_transform_iterator
                   (
                       f.begin(),
                       cmptMagUnaryFunctionFunctor<Type,Type>()
                   ),
                   thrust::make_transform_iterator
                   (
                       f.end(),
                       cmptMagUnaryFunctionFunctor<Type,Type>()
                   ),
                   pTraits<Type>::zero,
                   thrust::plus<Type>()
               );
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const gpuList<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }
    else
    {
        WarningIn("average(const gpuList<Type>&)")
            << "empty field, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                      \
                                                                              \
template<class Type>                                                          \
ReturnType gFunc(const gpuList<Type>& f, const int comm)                      \
{                                                                             \
    ReturnType res = Func(f);                                                 \
    reduce(res, rFunc##Op<Type>(), Pstream::msgType(), comm);                 \
    return res;                                                               \
}                                                                             \
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
    const int comm
)
{
    scalar SumProd = sumProd(f1, f2);
    reduce(SumProd, sumOp<scalar>(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type gSumCmptProd
(
    const gpuList<Type>& f1,
    const gpuList<Type>& f2,
    const int comm
)
{
    Type SumProd = sumCmptProd(f1, f2);
    reduce(SumProd, sumOp<Type>(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type gAverage
(
    const gpuList<Type>& f,
    const int comm
)
{
    label n = f.size();
    Type s = sum(f);
    sumReduce(s, n, Pstream::msgType(), comm);

    if (n > 0)
    {
        Type avrg = s/n;

        return avrg;
    }
    else
    {
        WarningIn("gAverage(const gpuList<Type>&)")
            << "empty field, returning zero." << endl;

        return pTraits<Type>::zero;
    }
}

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


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                    \
                                                                                 \
template<class Type1, class Type2>                                               \
void OpFunc                                                                      \
(                                                                                \
    gpuField<typename product<Type1, Type2>::type>& res,                         \
    const gpuList<Type1>& f1,                                                    \
    const gpuList<Type2>& f2                                                     \
)                                                                                \
{                                                                                \
    typedef typename product<Type1, Type2>::type productType;                    \
    thrust::transform(f1.begin(),f1.end(),f2.begin(),res.begin(),                \
            OpFunc##OperatorFunctor<Type1,Type2,productType>());                 \
}                                                                                \
                                                                                 \
template<class Type1, class Type2>                                               \
tmp<gpuField<typename product<Type1, Type2>::type> >                             \
operator Op(const gpuList<Type1>& f1, const gpuList<Type2>& f2)                  \
{                                                                                \
    typedef typename product<Type1, Type2>::type productType;                    \
    tmp<gpuField<productType> > tRes(new gpuField<productType>(f1.size()));      \
    OpFunc(tRes(), f1, f2);                                                      \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Type1, class Type2>                                               \
tmp<gpuField<typename product<Type1, Type2>::type> >                             \
operator Op(const gpuList<Type1>& f1, const tmp<gpuField<Type2> >& tf2)          \
{                                                                                \
    typedef typename product<Type1, Type2>::type productType;                    \
    tmp<gpuField<productType> > tRes = reusegpuTmp<productType, Type2>::New(tf2);\
    OpFunc(tRes(), f1, tf2());                                                   \
    reusegpuTmp<productType, Type2>::clear(tf2);                                 \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Type1, class Type2>                                               \
tmp<gpuField<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<gpuField<Type1> >& tf1, const gpuList<Type2>& f2)          \
{                                                                                \
    typedef typename product<Type1, Type2>::type productType;                    \
    tmp<gpuField<productType> > tRes = reusegpuTmp<productType, Type1>::New(tf1);\
    OpFunc(tRes(), tf1(), f2);                                                   \
    reusegpuTmp<productType, Type1>::clear(tf1);                                 \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Type1, class Type2>                                               \
tmp<gpuField<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<gpuField<Type1> >& tf1, const tmp<gpuField<Type2> >& tf2)  \
{                                                                                \
    typedef typename product<Type1, Type2>::type productType;                    \
    tmp<gpuField<productType> > tRes =                                           \
        reusegpuTmpTmp<productType, Type1, Type1, Type2>::New(tf1, tf2);         \
    OpFunc(tRes(), tf1(), tf2());                                                \
    reusegpuTmpTmp<productType, Type1, Type1, Type2>::clear(tf1, tf2);           \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Type, class Form, class Cmpt, int nCmpt>                          \
void OpFunc                                                                      \
(                                                                                \
    gpuField<typename product<Type, Form>::type>& res,                           \
    const gpuList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                                \
{                                                                                \
    typedef typename product<Type, Form>::type productType;                      \
    const Form& ff = static_cast<const Form&>(vs);                               \
    thrust::transform                                                            \
    (                                                                            \
        f1.begin(),                                                              \
        f1.end(),                                                                \
        res.begin(),                                                             \
        OpFunc##OperatorFSFunctor<Type,Form,productType>(ff)                     \
    );                                                                           \
}                                                                                \
                                                                                 \
template<class Type, class Form, class Cmpt, int nCmpt>                          \
tmp<gpuField<typename product<Type, Form>::type> >                               \
operator Op(const gpuList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                                \
    typedef typename product<Type, Form>::type productType;                      \
    tmp<gpuField<productType> > tRes(new gpuField<productType>(f1.size()));      \
    OpFunc(tRes(), f1, static_cast<const Form&>(vs));                            \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Type, class Form, class Cmpt, int nCmpt>                          \
tmp<gpuField<typename product<Type, Form>::type> >                               \
operator Op                                                                      \
(                                                                                \
    const tmp<gpuField<Type> >& tf1,                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                       \
)                                                                                \
{                                                                                \
    typedef typename product<Type, Form>::type productType;                      \
    tmp<gpuField<productType> > tRes = reusegpuTmp<productType, Type>::New(tf1); \
    OpFunc(tRes(), tf1(), static_cast<const Form&>(vs));                         \
    reusegpuTmp<productType, Type>::clear(tf1);                                  \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Form, class Cmpt, int nCmpt, class Type>                          \
void OpFunc                                                                      \
(                                                                                \
    gpuField<typename product<Form, Type>::type>& res,                           \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                      \
    const gpuList<Type>& f1                                                      \
)                                                                                \
{                                                                                \
    typedef typename product<Form, Type>::type productType;                      \
    const Form& ff = static_cast<const Form&>(vs);                               \
    thrust::transform                                                            \
    (                                                                            \
        f1.begin(),                                                              \
        f1.end(),                                                                \
        res.begin(),                                                             \
        OpFunc##OperatorSFFunctor<Form,Type,productType>(ff)                     \
    );                                                                           \
}                                                                                \
                                                                                 \
template<class Form, class Cmpt, int nCmpt, class Type>                          \
tmp<gpuField<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const gpuList<Type>& f1)     \
{                                                                                \
    typedef typename product<Form, Type>::type productType;                      \
    tmp<gpuField<productType> > tRes(new gpuField<productType>(f1.size()));      \
    OpFunc(tRes(), static_cast<const Form&>(vs), f1);                            \
    return tRes;                                                                 \
}                                                                                \
                                                                                 \
template<class Form, class Cmpt, int nCmpt, class Type>                          \
tmp<gpuField<typename product<Form, Type>::type> >                               \
operator Op                                                                      \
(                                                                                \
    const VectorSpace<Form,Cmpt,nCmpt>& vs, const tmp<gpuField<Type> >& tf1      \
)                                                                                \
{                                                                                \
    typedef typename product<Form, Type>::type productType;                      \
    tmp<gpuField<productType> > tRes = reusegpuTmp<productType, Type>::New(tf1); \
    OpFunc(tRes(), static_cast<const Form&>(vs), tf1());                         \
    reusegpuTmp<productType, Type>::clear(tf1);                                  \
    return tRes;                                                                 \
}

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

// ************************************************************************* //
