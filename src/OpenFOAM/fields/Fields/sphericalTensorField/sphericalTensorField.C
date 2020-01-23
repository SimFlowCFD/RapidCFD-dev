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

#include "sphericalTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tensorField& tf
)
{
    return sph(tf);
}

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<sphericalTensor> > ret =
        transformFieldMask<sphericalTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const symmTensorField& stf
)
{
    return sph(stf);
}

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<sphericalTensor> > ret =
        transformFieldMask<sphericalTensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#include "gpuFieldCommonFunctions.C"

#define TEMPLATE
#include "gpuFieldFunctionsM.C"
#include "gpuList.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template class gpuList<sphericalTensor>;
template class gpuField<sphericalTensor>;

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)

BINARY_SYM_OPERATOR(sphericalTensor, scalar, sphericalTensor, *, outer)
BINARY_SYM_FUNCTION(sphericalTensor, scalar, sphericalTensor, multiply)
BINARY_SYM_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)

BINARY_FULL_OPERATOR(sphericalTensor, sphericalTensor, sphericalTensor, +, add)
BINARY_FULL_OPERATOR(sphericalTensor, sphericalTensor, sphericalTensor, -, subtract)
BINARY_FULL_OPERATOR(sphericalTensor, sphericalTensor, sphericalTensor, &, dot)
BINARY_SYM_OPERATOR(vector, vector, sphericalTensor, &, dot)
BINARY_FULL_OPERATOR(scalar, sphericalTensor, sphericalTensor, &&, dotdot)



template<>
tmp<gpuField<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tensorgpuField& tf
)
{
    return sph(tf);
}

template<>
tmp<gpuField<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<tensorgpuField>& ttf
)
{
    tmp<gpuField<sphericalTensor> > ret =
        transformFieldMask<sphericalTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<gpuField<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const symmTensorgpuField& stf
)
{
    return sph(stf);
}

template<>
tmp<gpuField<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<symmTensorgpuField>& tstf
)
{
    tmp<gpuField<sphericalTensor> > ret =
        transformFieldMask<sphericalTensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"

// force instantiation
#define TEMPLATE template
#define FTYPE sphericalTensor
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"

