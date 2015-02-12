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

#include "nutUWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

struct nutUWallFunctionCalcNutFunctor
{
    const scalar yPlusLam;
    const scalar kappa;
    const scalar E;

    nutUWallFunctionCalcNutFunctor(const scalar yPlusLam_,const scalar kappa_,const scalar E_):
                                   yPlusLam(yPlusLam_),kappa(kappa_),E(E_){}

    __HOST____DEVICE__
    scalar operator () (const scalar& yPlus, const scalar& nuw)
    {	
        if (yPlus > yPlusLam)
        {
            return nuw*(yPlus*kappa/log(E*yPlus) - scalar(1.0));
        }
        else
        {
            return 0;
        }
    }
};

struct nutUWallFunctionCalcYPlusFunctor{
    const scalar yPlusLam;
    const scalar kappa;
    const scalar E;
    nutUWallFunctionCalcYPlusFunctor(const scalar yPlusLam_,const scalar kappa_,const scalar E_):
                                     yPlusLam(yPlusLam_),kappa(kappa_),E(E_){}

    __HOST____DEVICE__
    scalar operator () (const scalar& magUp, const thrust::tuple<scalar,scalar>& t)
    {
        scalar y = thrust::get<0>(t);
        scalar nuw = thrust::get<1>(t);
		
        scalar kappaRe = kappa*magUp*y/nuw;

        scalar yp = yPlusLam;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );

        return max(0.0, yp);
    }
};

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalargpuField> nutUWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalargpuField magUp(mag(Uw.patchInternalField() - Uw));
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalargpuField& nuw = nu.boundaryField()[patchi];

    tmp<scalargpuField> tyPlus = calcYPlus(magUp);
    scalargpuField& yPlus = tyPlus();

    tmp<scalargpuField> tnutw(new scalargpuField(patch().size(), 0.0));
    scalargpuField& nutw = tnutw();
    
    thrust::transform
    (
        yPlus.begin(),
        yPlus.end(),
        nuw.begin(),
        nutw.begin(),
        nutUWallFunctionCalcNutFunctor(yPlusLam_,kappa_,E_)
    );

/*
    forAll(yPlus, facei)
    {
        if (yPlus[facei] > yPlusLam_)
        {
            nutw[facei] =
                nuw[facei]*(yPlus[facei]*kappa_/log(E_*yPlus[facei]) - 1.0);
        }
    }
*/
    return tnutw;
}


tmp<scalargpuField> nutUWallFunctionFvPatchScalarField::calcYPlus
(
    const scalargpuField& magUp
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalargpuField& nuw = nu.boundaryField()[patchi];

    tmp<scalargpuField> tyPlus(new scalargpuField(patch().size(), 0.0));
    scalargpuField& yPlus = tyPlus();

    thrust::transform
    (
        magUp.begin(),
        magUp.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            y.begin(),
            nuw.begin()
        )),
        yPlus.begin(),
        nutUWallFunctionCalcYPlusFunctor(yPlusLam_,kappa_,E_)
    );

    /*
    forAll(yPlus, facei)
    {
        scalar kappaRe = kappa_*magUp[facei]*y[facei]/nuw[facei];

        scalar yp = yPlusLam_;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );

        yPlus[facei] = max(0.0, yp);
    }
    */
    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf)
{}


nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalargpuField> nutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalargpuField magUp(mag(Uw.patchInternalField() - Uw));

    return calcYPlus(magUp);
}


void nutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
