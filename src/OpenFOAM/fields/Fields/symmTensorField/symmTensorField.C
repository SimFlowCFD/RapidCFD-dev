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

#include "symmTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)

void inv(Field<symmTensor>& tf, const UList<symmTensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    scalar scale = magSqr(tf1[0]);
    Vector<bool> removeCmpts
    (
        magSqr(tf1[0].xx())/scale < SMALL,
        magSqr(tf1[0].yy())/scale < SMALL,
        magSqr(tf1[0].zz())/scale < SMALL
    );

    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        symmTensorField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += symmTensor(0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1Plus)

        if (removeCmpts.x())
        {
            tf -= symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= symmTensor(0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1)
    }
}

tmp<symmTensorField> inv(const UList<symmTensor>& tf)
{
    tmp<symmTensorField> result(new symmTensorField(tf.size()));
    inv(result(), tf);
    return result;
}

tmp<symmTensorField> inv(const tmp<symmTensorField>& tf)
{
    tmp<symmTensorField> tRes = reuseTmp<symmTensor, symmTensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<symmTensor, symmTensor>::clear(tf);
    return tRes;
}


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tensorField& tf
)
{
    return symm(tf);
}

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<symmTensor> > ret = transformFieldMask<symmTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const symmTensorField& stf
)
{
    return stf;
}

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    return tstf;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"


#define TEMPLATE
#include "gpuFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template class gpuField<symmTensor>;

struct symmTensorRemoveComponentsFunctor 
    : public thrust::unary_function<symmTensor, Vector<bool> >
{
    const scalar scale;
    symmTensorRemoveComponentsFunctor(scalar _scale): scale(_scale) {}
    __HOST____DEVICE__
    Vector<bool> operator()(const symmTensor& st) const 
    {
        return Vector<bool>
        (
            st.xx()/scale < SMALL,
            st.yy()/scale < SMALL,
            st.zz()/scale < SMALL
        );
    }
};

struct andBooleanVectorFunctor 
    : public thrust::binary_function<Vector<bool>,Vector<bool>,Vector<bool> >
{
    __HOST____DEVICE__
    Vector<bool> operator()(const Vector<bool>& v1,const Vector<bool>& v2) const 
    {
        return Vector<bool>
        (
            v1.x()&&v2.x(),
            v1.y()&&v2.y(),
            v1.z()&&v2.z()
        );
    }
};

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)

void inv(gpuField<symmTensor>& tf, const gpuList<symmTensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    gpuList<symmTensor> tmp(tf1,1);
    scalar scale = sum(magSqr(tmp));
    Vector<bool> removeCmpts = thrust::reduce
        (
            thrust::make_transform_iterator
            (
                tmp.begin(),
                symmTensorRemoveComponentsFunctor(scale)
            ),
            thrust::make_transform_iterator
            (
                tmp.end(),
                symmTensorRemoveComponentsFunctor(scale)
            ),
            Vector<bool>(true,true,true),
            andBooleanVectorFunctor()
        );

    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        symmTensorgpuField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += symmTensor(0,0,0,0,0,1);
        }

        thrust::transform
        (
            tf1Plus.begin(),
            tf1Plus.end(),
            tf.begin(),
            invUnaryFunctionFunctor<symmTensor,symmTensor>()
        );

        if (removeCmpts.x())
        {
            tf -= symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= symmTensor(0,0,0,0,0,1);
        }
    }
    else
    {
        thrust::transform
        (
            tf1.begin(),
            tf1.end(),
            tf.begin(),
            invUnaryFunctionFunctor<symmTensor,symmTensor>()
        );
    }
}

tmp<symmTensorgpuField> inv(const gpuList<symmTensor>& tf)
{
    tmp<symmTensorgpuField> result(new symmTensorgpuField(tf.size()));
    inv(result(), tf);
    return result;
}

tmp<symmTensorgpuField> inv(const tmp<symmTensorgpuField>& tf)
{
    tmp<symmTensorgpuField> tRes = reuseTmp<symmTensor, symmTensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<symmTensor, symmTensor>::clear(tf);
    return tRes;
}


template<>
tmp<gpuField<symmTensor> > transformFieldMask<symmTensor>
(
    const tensorgpuField& tf
)
{
    return symm(tf);
}

template<>
tmp<gpuField<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<tensorgpuField>& ttf
)
{
    tmp<gpuField<symmTensor> > ret = transformFieldMask<symmTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<gpuField<symmTensor> > transformFieldMask<symmTensor>
(
    const symmTensorgpuField& stf
)
{
    return stf;
}

template<>
tmp<gpuField<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<symmTensorgpuField>& tstf
)
{
    return tstf;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"

// ************************************************************************* //
