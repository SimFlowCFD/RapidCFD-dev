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

#include "gpuFieldCommonFunctionsM.H"

namespace Foam
{

template<class Type, int r>
void pow
(
    gpuField<typename powProduct<Type, r>::type>& res,
    const gpuList<Type>& vf
);


template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const gpuList<Type>& f,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, int r>
tmp<gpuField<typename powProduct<Type, r>::type> >
pow
(
    const tmp<gpuField<Type> >& tf,
    typename powProduct<Type, r>::type
      = pTraits<typename powProduct<Type, r>::type>::zero
);

}

#include "scalarField.H"

