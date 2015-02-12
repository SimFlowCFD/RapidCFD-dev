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

#include "nutkAtmRoughWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
	
struct nutkAtmRoughWallFunctionCalcNutFunctor
{
    const scalar Cmu25;
    const scalar kappa;

    nutkAtmRoughWallFunctionCalcNutFunctor(const scalar Cmu25_,const scalar kappa_):Cmu25(Cmu25_),kappa(kappa_){}

    __HOST____DEVICE__
    scalar operator () (const scalar& y, const thrust::tuple<scalar,scalar,scalar>& t)
    {
        scalar k = thrust::get<0>(t);
        scalar nuw = thrust::get<1>(t);
        scalar z0 = thrust::get<2>(t);

        scalar uStar = Cmu25*sqrt(k);
        scalar yPlus = uStar*y/nuw;

        scalar Edash = (y + z0)/(z0 + 1e-4);

        return nuw*(yPlus*kappa/log(max(Edash, 1 + 1e-4)) - 1);
    }
};

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalargpuField> nutkAtmRoughWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbulence.y()[patchI];
    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> tnu = turbulence.nu();
    const volScalarField& nu = tnu();
    const scalargpuField& nuw = nu.boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalargpuField> tnutw(new scalargpuField(*this));
    scalargpuField& nutw = tnutw();
    
    thrust::transform
    (
        y.begin(),
        y.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator(k.getField().begin(),patch().faceCells().begin()),
            nuw.begin(),
            z0_.begin()
        )),
        nutw.begin(),
        nutkAtmRoughWallFunctionCalcNutFunctor(Cmu25,kappa_)
    );

/*
    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/nuw[faceI];

        scalar Edash = (y[faceI] + z0_[faceI])/(z0_[faceI] + 1e-4);

        nutw[faceI] =
            nuw[faceI]*(yPlus*kappa_/log(max(Edash, 1 + 1e-4)) - 1);

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[faceI]
                << endl;
        }
    }
*/
    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    z0_(p.size(), 0.0)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_, mapper)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    z0_("z0", dict, p.size())
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    z0_(rwfpsf.z0_)
{}


nutkAtmRoughWallFunctionFvPatchScalarField::
nutkAtmRoughWallFunctionFvPatchScalarField
(
    const nutkAtmRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    z0_(rwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutkAtmRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
}


void nutkAtmRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelgpuList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkAtmRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkAtmRoughWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}


void nutkAtmRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    z0_.writeEntry("z0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkAtmRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
