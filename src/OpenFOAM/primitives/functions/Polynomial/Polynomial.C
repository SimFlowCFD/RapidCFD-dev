/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const UList<scalar>& coeffs)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    if (coeffs.size() != PolySize)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const UList<scalar>&)"
        )   << "Size mismatch: Needed " << PolySize
            << " but given " << coeffs.size()
            << nl << exit(FatalError);
    }

    for (int i = 0; i < PolySize; ++i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is),
    logActive_(false),
    logCoeff_(0.0)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const word& name, Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    word isName(is);

    if (isName != name)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Expected polynomial name " << name << " but read " << isName
            << nl << exit(FatalError);
    }

    VectorSpace<Polynomial<PolySize>, scalar, PolySize>::
        operator=(VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is));

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Polynomial coefficients for entry " << isName
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


// ************************************************************************* //
