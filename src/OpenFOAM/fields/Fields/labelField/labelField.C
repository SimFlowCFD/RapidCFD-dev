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
    Specialisation of Field\<T\> for label.

\*---------------------------------------------------------------------------*/

#include "labelField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)


template<>
tmp<labelField> labelField::component(const direction) const
{
    return *this;
}

template<>
void component
(
    labelField& lf,
    const labelUList& f,
    const direction
)
{
    lf = f;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

#include "gpuFieldCommonFunctions.C"
namespace Foam
{

template
label max<label>(const gpuList<label>&);


template<>
void labelField::replace(const direction, const labelUList& lf)
{
    *this = lf;
}

template<>
__host__ __device__
label transposeFunctor<label>::operator()(const label& s) const
{
    return s;
}


template<>
struct componentFunctor<label,label>
{
    const direction d;
    componentFunctor(direction _d): d(_d) {}
    __host__ __device__
    label operator()(const label& tt)
    {
        return tt;
    }
};
}

// force instantiation
//#define TEMPLATE template
//#define FTYPE label
//#include "gpuFieldCommonFunctionsM.H"

// ************************************************************************* //
