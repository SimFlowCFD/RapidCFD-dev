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

#include "polynomial.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "textures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polynomial, 0);
    addToRunTimeSelectionTable(scalarDataEntry, polynomial, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polynomial::polynomial(const word& entryName, const dictionary& dict)
:
    scalarDataEntry(entryName),
    coeffs_(),
    preCoeffs_(coeffs_.size()),
    expCoeffs_(coeffs_.size()),
    canIntegrate_(true),
    dimensions_(dimless)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);

    token firstToken(is);
    is.putBack(firstToken);
    if (firstToken == token::BEGIN_SQR)
    {
        is  >> this->dimensions_;
    }

    is  >> coeffs_;

    if (!coeffs_.size())
    {
        FatalErrorIn
        (
            "Foam::polynomial::polynomial(const word&, const dictionary&)"
        )   << "polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + 1) < ROOTVSMALL)
        {
            canIntegrate_ = false;
            break;
        }
    }

    if (debug)
    {
        if (!canIntegrate_)
        {
            WarningIn
            (
                "Foam::polynomial::polynomial(const word&, const dictionary&)"
            )   << "Polynomial " << this->name_ << " cannot be integrated"
                << endl;
        }
    }

    initGpuCoeffs();
}


Foam::polynomial::polynomial
(
    const word& entryName,
    const List<Tuple2<scalar, scalar> >& coeffs
)
:
    scalarDataEntry(entryName),
    coeffs_(coeffs),
    preCoeffs_(coeffs_.size()),
    expCoeffs_(coeffs_.size()),
    canIntegrate_(true),
    dimensions_(dimless)
{
    if (!coeffs_.size())
    {
        FatalErrorIn
        (
            "Foam::polynomial::polynomial"
            "(const word&, const List<Tuple2<scalar, scalar> >&)"
        )   << "polynomial coefficients for entry " << this->name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }

    forAll(coeffs_, i)
    {
        if (mag(coeffs_[i].second() + 1) < ROOTVSMALL)
        {
            canIntegrate_ = false;
            break;
        }
    }

    if (debug)
    {
        if (!canIntegrate_)
        {
            WarningIn
            (
                "Foam::polynomial::polynomial"
                "(const word&, const List<Tuple2<scalar, scalar> >&)"
            )   << "Polynomial " << this->name_ << " cannot be integrated"
                << endl;
        }
    }

    initGpuCoeffs();
}


Foam::polynomial::polynomial(const polynomial& poly)
:
    scalarDataEntry(poly),
    coeffs_(poly.coeffs_),
    preCoeffs_(coeffs_.size()),
    expCoeffs_(coeffs_.size()),
    canIntegrate_(poly.canIntegrate_),
    dimensions_(poly.dimensions_)
{
    initGpuCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polynomial::~polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polynomial::initGpuCoeffs()
{
    scalarList pre(coeffs_.size());
    scalarList exp(coeffs_.size());

    forAll(coeffs_, i)
    {
        pre[i] = coeffs_[i].first();
        exp[i] = coeffs_[i].second();
    }

    preCoeffs_ = pre;
    expCoeffs_ = exp;
}


void Foam::polynomial::convertTimeBase(const Time& t)
{
    forAll(coeffs_, i)
    {
        scalar value = coeffs_[i].first();
        coeffs_[i].first() = t.userTimeToTime(value);
    }
}


Foam::scalar Foam::polynomial::value(const scalar x) const
{
    scalar y = 0.0;
    forAll(coeffs_, i)
    {
        y += coeffs_[i].first()*pow(x, coeffs_[i].second());
    }

    return y;
}


Foam::scalar Foam::polynomial::integrate(const scalar x1, const scalar x2) const
{
    scalar intx = 0.0;

    if (canIntegrate_)
    {
        forAll(coeffs_, i)
        {
            intx +=
                coeffs_[i].first()/(coeffs_[i].second() + 1)
               *(
                    pow(x2, coeffs_[i].second() + 1)
                  - pow(x1, coeffs_[i].second() + 1)
                );
        }
    }

    return intx;
}


Foam::dimensioned<Foam::scalar> Foam::polynomial::dimValue
(
    const scalar x
) const
{
    return dimensioned<scalar>("dimensionedValue", dimensions_, value(x));
}


Foam::dimensioned<Foam::scalar> Foam::polynomial::dimIntegrate
(
    const scalar x1,
    const scalar x2
) const
{
    return dimensioned<scalar>
    (
        "dimensionedValue",
        dimensions_,
        integrate(x1, x2)
    );
}


namespace Foam
{

struct polynomialFunctor
{
    const label size_;
    const textures<scalar> pre_;
    const textures<scalar> exp_;

    polynomialFunctor
    (
        const label size,
        const textures<scalar> preCoeffs,
        const textures<scalar> expCoeffs
    ):
        size_(size),
        pre_(preCoeffs),
        exp_(expCoeffs)
    {}

    __device__
    scalar operator()(const scalar& x)
    {
        scalar y = 0.0;

        for(label i = 0; i < size_; i++)
        {
            y += pre_[i]*pow(x, exp_[i]);
        }

        return y;
    }
};

}


Foam::tmp<Foam::scalargpuField > Foam::polynomial::value
(
    const scalargpuField& x
) const
{
    tmp<scalargpuField> tfld(new scalargpuField(x.size()));
    scalargpuField& fld = tfld();

    textures<scalar> preCoeffs(preCoeffs_);
    textures<scalar> expCoeffs(expCoeffs_);

    thrust::transform
    (
        x.begin(),
        x.end(),
        fld.begin(),
        polynomialFunctor
        (
            coeffs_.size(),
            preCoeffs,
            expCoeffs
        )
    );

    preCoeffs.destroy();
    expCoeffs.destroy();

    return tfld;
}


// ************************************************************************* //
