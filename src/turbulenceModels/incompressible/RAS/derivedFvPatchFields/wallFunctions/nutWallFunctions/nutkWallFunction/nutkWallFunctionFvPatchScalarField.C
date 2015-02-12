/*---------------------------------------------------------------------------* \
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

#include "nutkWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
	
struct nutkWallFunctionFvPatchScalarFieldFunctor
{
    const scalar yPlusLam_;
    const scalar Cmu25_;
    const scalar kappa_;
    const scalar E_;

    nutkWallFunctionFvPatchScalarFieldFunctor
    (
        scalar yPlusLam,
        scalar Cmu25,
        scalar kappa,
        scalar E
    ):
	yPlusLam_(yPlusLam),
	Cmu25_(Cmu25),
	kappa_(kappa),
	E_(E)
    {}

    __HOST____DEVICE__
    scalar operator () (const thrust::tuple<scalar,scalar,scalar>& t)
    {
        scalar y = thrust::get<0>(t);
        scalar k = thrust::get<1>(t);
        scalar nuw = thrust::get<2>(t);
		
        scalar yPlus = Cmu25_*y*sqrt(k)/nuw;

        if (yPlus > yPlusLam_)
        {
            return nuw*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
        else
        {
            return 0;
        }
    }	
};

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalargpuField> nutkWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalargpuField& nuw = nu.boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalargpuField> tnutw(new scalargpuField(patch().size()));
    scalargpuField& nutw = tnutw();

    thrust::transform
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            y.begin(),
            thrust::make_permutation_iterator
            (
                k.getField().begin(),
                patch().faceCells().begin()
            ),
            nuw.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            y.end(),
            thrust::make_permutation_iterator
            (
                k.getField().begin(),
                patch().faceCells().end()
            ),
            nuw.end()
        )),
        nutw.begin(),
        nutkWallFunctionFvPatchScalarFieldFunctor(yPlusLam_,Cmu25,kappa_,E_)
    );

/*    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus = Cmu25*y[faceI]*sqrt(k[faceCellI])/nuw[faceI];

        if (yPlus > yPlusLam_)
        {
            nutw[faceI] = nuw[faceI]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
    }*/

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalargpuField> nutkWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalargpuField> kwc = k.boundaryField()[patchi].patchInternalField();
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalargpuField& nuw = nu.boundaryField()[patchi];

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutkWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
