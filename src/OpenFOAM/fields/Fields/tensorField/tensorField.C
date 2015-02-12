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

#include "tensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)

void inv(Field<tensor>& tf, const UList<tensor>& tf1)
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
        tensorField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1Plus)

        if (removeCmpts.x())
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1)
    }
}

tmp<tensorField> inv(const UList<tensor>& tf)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result(), tf);
    return result;
}

tmp<tensorField> inv(const tmp<tensorField>& tf)
{
    tmp<tensorField> tRes = reuseTmp<tensor, tensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<tensor, tensor>::clear(tf);
    return tRes;
}

UNARY_FUNCTION(vector, tensor, eigenValues)
UNARY_FUNCTION(tensor, tensor, eigenVectors)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const symmTensorField& stf
)
{
    tmp<tensorField> tRes(new tensorField(stf.size()));
    tensorField& res = tRes();
    TFOR_ALL_F_OP_F(tensor, res, =, symmTensor, stf)
    return tRes;
}

template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<tensor> > ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //

#define TEMPLATE
#include "gpuFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template class gpuField<tensor>;

struct tensorRemoveComponentsFunctor : public thrust::unary_function<symmTensor, Vector<bool> >
{
    const scalar scale;
    tensorRemoveComponentsFunctor(scalar _scale): scale(_scale) {}
    __HOST____DEVICE__
    Vector<bool> operator()(const tensor& st) const 
    {
        return Vector<bool>
        (
            st.xx()/scale < SMALL,
            st.yy()/scale < SMALL,
            st.zz()/scale < SMALL
        );
    }
};

struct andBooleanVectorTensorFunctor : public thrust::binary_function<Vector<bool>,Vector<bool>,Vector<bool> >
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

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)

void inv(gpuField<tensor>& tf, const gpuList<tensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    gpuList<tensor> tmp(tf1,1);
    scalar scale = sum(magSqr(tmp));

    Vector<bool> removeCmpts = 
        thrust::reduce
        (
            thrust::make_transform_iterator
            (
                tmp.begin(),
                tensorRemoveComponentsFunctor(scale)
            ),
            thrust::make_transform_iterator
            (
                tmp.end(),
                tensorRemoveComponentsFunctor(scale)
            ),
            Vector<bool>(true,true,true),
            andBooleanVectorTensorFunctor()
        );

    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        tensorgpuField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        thrust::transform
        (
            tf1Plus.begin(),
            tf1Plus.end(),
            tf.begin(),
            invUnaryFunctionFunctor<tensor,tensor>()
        );

        if (removeCmpts.x())
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        thrust::transform
        (
            tf1.begin(),
            tf1.end(),
            tf.begin(),
            invUnaryFunctionFunctor<tensor,tensor>()
        );
    }
}

tmp<tensorgpuField> inv(const gpuList<tensor>& tf)
{
    tmp<tensorgpuField> result(new tensorgpuField(tf.size()));
    inv(result(), tf);
    return result;
}

tmp<tensorgpuField> inv(const tmp<tensorgpuField>& tf)
{
    tmp<tensorgpuField> tRes = reuseTmp<tensor, tensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<tensor, tensor>::clear(tf);
    return tRes;
}

UNARY_FUNCTION(vector, tensor, eigenValues)
UNARY_FUNCTION(tensor, tensor, eigenVectors)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)

template<>
tmp<gpuField<tensor> > transformFieldMask<tensor>
(
    const symmTensorgpuField& stf
)
{
    tmp<tensorgpuField> tRes(new tensorgpuField(stf.size()));
    tensorgpuField& res = tRes();
    thrust::transform
    (
        stf.begin(),
        stf.end(),
        res.begin(),
        assignFunctor<symmTensor,tensor>()
    );
    return tRes;
}

template<>
tmp<gpuField<tensor> > transformFieldMask<tensor>
(
    const tmp<symmTensorgpuField>& tstf
)
{
    tmp<gpuField<tensor> > ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"
