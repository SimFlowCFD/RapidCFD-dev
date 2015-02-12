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

#include "mutUWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

struct mutUCalcMutFunctor
{
    const scalar yPlusLam;
    const scalar kappa;
    const scalar E;

    mutUCalcMutFunctor(const scalar yPlusLam_,const scalar kappa_,const scalar E_):
                       yPlusLam(yPlusLam_),kappa(kappa_),E(E_){}

    __HOST____DEVICE__
    scalar operator () (const scalar& yPlus, const scalar& muw)
    {	
        if (yPlus > yPlusLam)
        {
            return muw*(yPlus*kappa/log(E*yPlus) - 1.0);
        }
        else
        {
            return 0;
        }
    }
};

struct mutUCalcYPlusFunctor
{
    const scalar yPlusLam;
    const scalar kappa;
    const scalar E;

    mutUCalcYPlusFunctor(const scalar yPlusLam_,const scalar kappa_,const scalar E_):
                         yPlusLam(yPlusLam_),kappa(kappa_),E(E_){}
    __HOST____DEVICE__
    scalar operator () (const scalar& magUp, const thrust::tuple<scalar,scalar,scalar>& t)
    {
        scalar y = thrust::get<0>(t);
        scalar muw = thrust::get<1>(t);
        scalar rhow = thrust::get<2>(t);
		
        scalar kappaRe = kappa*magUp*y/(muw/rhow);

        scalar yp = yPlusLam;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10);

        return max(0.0, yp);
	}
};

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalargpuField> mutUWallFunctionFvPatchScalarField::calcYPlus
(
    const scalargpuField& magUp
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalargpuField& y = turbModel.y()[patchi];
    const fvPatchScalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& muw = turbModel.mu().boundaryField()[patchi];

    tmp<scalargpuField> tyPlus(new scalargpuField(patch().size(), 0.0));
    scalargpuField& yPlus = tyPlus();
    
    thrust::transform
    (
        magUp.begin(),
        magUp.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            y.begin(),
            muw.begin(),
            rhow.begin()
        )),
        yPlus.begin(),
        mutUCalcYPlusFunctor(yPlusLam_,kappa_,E_)
    );
/*
    forAll(yPlus, faceI)
    {
        scalar kappaRe = kappa_*magUp[faceI]*y[faceI]/(muw[faceI]/rhow[faceI]);

        scalar yp = yPlusLam_;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10);

        yPlus[faceI] = max(0.0, yp);
    }
*/
    return tyPlus;
}


tmp<scalargpuField> mutUWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalargpuField magUp(mag(Uw.patchInternalField() - Uw));
    const fvPatchScalarField& muw = turbModel.mu().boundaryField()[patchi];

    tmp<scalargpuField> tyPlus = calcYPlus(magUp);
    scalargpuField& yPlus = tyPlus();

    tmp<scalargpuField> tmutw(new scalargpuField(patch().size(), 0.0));
    scalargpuField& mutw = tmutw();
    
    thrust::transform
    (
        yPlus.begin(),
        yPlus.end(),
        muw.begin(),
        mutw.begin(),
        mutUCalcMutFunctor(yPlusLam_,kappa_,E_)
    );

/*
    forAll(yPlus, faceI)
    {
        if (yPlus[faceI] > yPlusLam_)
        {
            mutw[faceI] =
                muw[faceI]*(yPlus[faceI]*kappa_/log(E_*yPlus[faceI]) - 1.0);
        }
    }
*/
    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    mutWallFunctionFvPatchScalarField(sawfpsf)
{}


mutUWallFunctionFvPatchScalarField::mutUWallFunctionFvPatchScalarField
(
    const mutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalargpuField> mutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalargpuField magUp(mag(Uw.patchInternalField() - Uw));

    return calcYPlus(magUp);
}


void mutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
