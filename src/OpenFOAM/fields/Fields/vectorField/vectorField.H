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
    Foam::vectorField

Description
    Specialisation of Field\<T\> for vector.

SourceFiles
    vectorField.C

\*---------------------------------------------------------------------------*/

#ifndef vectorField_H
#define vectorField_H

#include "scalarField.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TEMPLATE
#include "gpuFieldFunctionsM.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<vector> vectorField;
typedef gpuField<vector> vectorgpuField;

BINARY_SYM_OPERATOR(vector, scalar, vector, *, outer)
BINARY_SYM_FUNCTION(vector, scalar, vector, multiply)
BINARY_OPERATOR(vector, vector, scalar, /, divide)
BINARY_TYPE_OPERATOR_FS(vector, vector, scalar, /, divide)

BINARY_FULL_OPERATOR(vector, vector, vector, +, add)
BINARY_FULL_OPERATOR(vector, vector, vector, -, subtract)

BINARY_FULL_OPERATOR(scalar, vector, vector, &, dot)
BINARY_FULL_OPERATOR(vector, vector, vector, ^, cross)

extern template class gpuList<vector>;
extern template class gpuField<vector>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "undefgpuFieldFunctionsM.H"

#define TEMPLATE extern template
#define FTYPE vector
#define NO_SQR
#include "gpuFieldCommonFunctionsM.H"

#endif

// ************************************************************************* //
