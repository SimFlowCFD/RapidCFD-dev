
#include "UPstream.H"

#ifdef TEMPLATE
  #ifndef FTYPE
    #error "FTYPE macro is not defined"
  #endif
  #define PSTREAM_DEFAULT
#else
  #define TEMPLATE template<class Type>
  #define FTYPE Type
  #define PSTREAM_DEFAULT = UPstream::worldComm
#endif

#include "gpuFieldFunctionsM.H"

namespace Foam
{

TEMPLATE
void component
(
    gpuField<typename gpuField<FTYPE>::cmptType>& res,
    const gpuList<FTYPE>& f,
    const direction d
);


TEMPLATE
void T(gpuField<FTYPE>& res, const gpuList<FTYPE>& f);


TEMPLATE
void magSqr(gpuField<scalar>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<scalar> > magSqr(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<scalar> > magSqr(const tmp<gpuField<FTYPE> >& tf);


TEMPLATE
void mag(gpuField<scalar>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<scalar> > mag(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<scalar> > mag(const tmp<gpuField<FTYPE> >& tf);


TEMPLATE
void cmptMax(gpuField<typename gpuField<FTYPE>::cmptType>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> > cmptMax(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> >
cmptMax(const tmp<gpuField<FTYPE> >& tf);


TEMPLATE
void cmptMin(gpuField<typename gpuField<FTYPE>::cmptType>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> > cmptMin(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> >
cmptMin(const tmp<gpuField<FTYPE> >& tf);


TEMPLATE
void cmptAv(gpuField<typename gpuField<FTYPE>::cmptType>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> > cmptAv(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<typename gpuField<FTYPE>::cmptType> > cmptAv(const tmp<gpuField<FTYPE> >& tf);


TEMPLATE
void cmptMag(gpuField<FTYPE>& res, const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<FTYPE> > cmptMag(const gpuList<FTYPE>& f);

TEMPLATE
tmp<gpuField<FTYPE> > cmptMag(const tmp<gpuField<FTYPE> >& tf);

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
                                                                              \
TEMPLATE                                                                      \
ReturnType Func(const tmp<gpuField<FTYPE> >& tf1);

TEMPLATE
FTYPE max(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, max)

TEMPLATE
FTYPE min(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, min)

TEMPLATE
FTYPE sum(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, sum)

TEMPLATE
FTYPE maxMagSqr(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, maxMagSqr)

TEMPLATE
FTYPE minMagSqr(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, minMagSqr)


TEMPLATE
scalar sumProd(const gpuList<FTYPE>& f1, const gpuList<FTYPE>& f2);

TEMPLATE
FTYPE sumCmptProd(const gpuList<FTYPE>& f1, const gpuList<FTYPE>& f2);

#ifndef NO_SQR
TEMPLATE
scalar sumSqr(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(scalar, sumSqr)
#endif

TEMPLATE
scalar sumMag(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

TEMPLATE
FTYPE sumCmptMag(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, sumCmptMag)

TEMPLATE
FTYPE average(const gpuList<FTYPE>& f);

TMP_UNARY_FUNCTION(FTYPE, average)

#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                        \
                                                                                \
TEMPLATE                                                                        \
ReturnType gFunc(const gpuList<FTYPE>& f, const int comm PSTREAM_DEFAULT); \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(FTYPE, gMax, max, max)
G_UNARY_FUNCTION(FTYPE, gMin, min, min)
G_UNARY_FUNCTION(FTYPE, gSum, sum, sum)
G_UNARY_FUNCTION(FTYPE, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(FTYPE, gMinMagSqr, minMagSqr, minMagSqr)
#ifndef NO_SQR
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
#endif
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(FTYPE, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

TEMPLATE
scalar gSumProd
(
    const gpuList<FTYPE>& f1,
    const gpuList<FTYPE>& f2,
    const int comm PSTREAM_DEFAULT
);

TEMPLATE
FTYPE gSumCmptProd
(
    const gpuList<FTYPE>& f1,
    const gpuList<FTYPE>& f2,
    const int comm PSTREAM_DEFAULT
);

TEMPLATE
FTYPE gAverage
(
    const gpuList<FTYPE>& f,
    const int comm PSTREAM_DEFAULT
);

TMP_UNARY_FUNCTION(FTYPE, gAverage)

#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(FTYPE, FTYPE, FTYPE, max)
BINARY_FUNCTION(FTYPE, FTYPE, FTYPE, min)
BINARY_FUNCTION(FTYPE, FTYPE, FTYPE, cmptMultiply)
BINARY_FUNCTION(FTYPE, FTYPE, FTYPE, cmptDivide)

BINARY_TYPE_FUNCTION(FTYPE, FTYPE, FTYPE, max)
BINARY_TYPE_FUNCTION(FTYPE, FTYPE, FTYPE, min)
BINARY_TYPE_FUNCTION(FTYPE, FTYPE, FTYPE, cmptMultiply)
BINARY_TYPE_FUNCTION(FTYPE, FTYPE, FTYPE, cmptDivide)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(FTYPE, FTYPE, -, negate)

}

#include "undefgpuFieldFunctionsM.H"
