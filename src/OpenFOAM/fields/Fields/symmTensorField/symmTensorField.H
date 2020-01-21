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
    Foam::symmTensorField

Description
    Specialisation of Field\<T\> for symmTensor.

SourceFiles
    symmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensorField_H
#define symmTensorField_H

#include "scalarField.H"
#include "vectorField.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<symmTensor> symmTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


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

typedef gpuField<symmTensor> symmTensorgpuField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_SYM_OPERATOR(symmTensor, scalar, symmTensor, *, outer)
BINARY_SYM_FUNCTION(symmTensor, scalar, symmTensor, multiply)
BINARY_OPERATOR(symmTensor, symmTensor, scalar, /, divide)
BINARY_TYPE_OPERATOR_FS(symmTensor, symmTensor, scalar, /, divide)

BINARY_FULL_OPERATOR(symmTensor, symmTensor, symmTensor, +, add)
BINARY_FULL_OPERATOR(symmTensor, symmTensor, symmTensor, -, subtract)

BINARY_FULL_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_SYM_OPERATOR(vector, vector, symmTensor, &, dot)
BINARY_FULL_OPERATOR(scalar, symmTensor, symmTensor, &&, dotdot)

BINARY_SYM_OPERATOR(symmTensor, sphericalTensor, symmTensor, +, add)
BINARY_SYM_OPERATOR(symmTensor, sphericalTensor, symmTensor, -, subtract)
BINARY_SYM_OPERATOR(symmTensor, sphericalTensor, symmTensor, &, dot)
BINARY_SYM_OPERATOR(scalar, sphericalTensor, symmTensor, &&, dotdot)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern template class gpuList<symmTensor>;
extern template class gpuField<symmTensor>;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefgpuFieldFunctionsM.H"

#define TEMPLATE extern template
#define FTYPE symmTensor
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"

#endif

// ************************************************************************* //
