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

Typedef
    Foam::sphericalTensorField

Description
    Specialisation of Field\<T\> for sphericalTensor.

SourceFiles
    sphericalTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensorField_H
#define sphericalTensorField_H

#include "scalarField.H"
#include "sphericalTensor.H"
#include "vector.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<sphericalTensor> sphericalTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TEMPLATE
#include "gpuFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef gpuField<sphericalTensor> sphericalTensorgpuField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern template class gpuList<sphericalTensor>;
extern template class gpuField<sphericalTensor>;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"

#define TEMPLATE extern template
#define FTYPE sphericalTensor
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"

#endif

// ************************************************************************* //
