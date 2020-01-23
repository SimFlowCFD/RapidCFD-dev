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
    Specialisation of Field\<T\> for diagTensor.

\*---------------------------------------------------------------------------*/

#include "diagTensorField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(diagTensor, tensor, diag)
UNARY_FUNCTION(scalar, diagTensor, tr)
UNARY_FUNCTION(sphericalTensor, diagTensor, sph)
UNARY_FUNCTION(scalar, diagTensor, det)
UNARY_FUNCTION(diagTensor, diagTensor, inv)


BINARY_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_OPERATOR(vector, vector, diagTensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, diagTensor, /, divide)


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

template class gpuList<diagTensor>;
template class gpuField<diagTensor>;

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(diagTensor, tensor, diag)
UNARY_FUNCTION(scalar, diagTensor, tr)
UNARY_FUNCTION(sphericalTensor, diagTensor, sph)
UNARY_FUNCTION(scalar, diagTensor, det)
UNARY_FUNCTION(diagTensor, diagTensor, inv)

BINARY_SYM_OPERATOR(diagTensor, scalar, diagTensor, *, outer)
BINARY_SYM_FUNCTION(diagTensor, scalar, diagTensor, multiply)
BINARY_OPERATOR(diagTensor, diagTensor, scalar, /, divide)
BINARY_TYPE_OPERATOR_FS(diagTensor, diagTensor, scalar, /, divide)

BINARY_FULL_OPERATOR(diagTensor, diagTensor, diagTensor, +, add)
BINARY_FULL_OPERATOR(diagTensor, diagTensor, diagTensor, -, subtract)
BINARY_FULL_OPERATOR(diagTensor, diagTensor, diagTensor, &, dot)
BINARY_FULL_OPERATOR(vector, vector, diagTensor, /, divide)
BINARY_FULL_OPERATOR(diagTensor, scalar, diagTensor, /, divide)

BINARY_SYM_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_SYM_OPERATOR(tensor, diagTensor, tensor, -, subtract)
BINARY_SYM_OPERATOR(tensor, diagTensor, tensor, &, dot)
BINARY_SYM_OPERATOR(vector, diagTensor, vector, &, dot)

} 


#include "undefgpuFieldFunctionsM.H"

// force instantiation
#define TEMPLATE template
#define FTYPE diagTensor
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"

