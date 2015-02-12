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

#include "mutkWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

struct mutkCalcMutFunctor
{
    const scalar yPlusLam_;
    const scalar Cmu25_;
    const scalar kappa_;
    const scalar E_;

    mutkCalcMutFunctor(scalar yPlusLam,scalar Cmu25,scalar kappa,scalar E):
	               yPlusLam_(yPlusLam),Cmu25_(Cmu25),kappa_(kappa),E_(E){}

    __HOST____DEVICE__
    scalar operator () (const thrust::tuple<scalar,scalar,scalar,scalar>& t)
    {
        scalar y = thrust::get<0>(t);
        scalar k = thrust::get<1>(t);
        scalar muw = thrust::get<2>(t);
        scalar rhow = thrust::get<3>(t);
		
        scalar yPlus = Cmu25_*y*sqrt(k)/(muw/rhow);

        if (yPlus > yPlusLam_)
        {
            return muw*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
        else
        {
            return 0;
        }
    }	
};

tmp<scalargpuField> mutkWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbModel.y()[patchi];
    const scalargpuField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalargpuField& muw = turbModel.mu().boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalargpuField> tmutw(new scalargpuField(patch().size(), 0.0));
    scalargpuField& mutw = tmutw();

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
            muw.begin(),
            rhow.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            y.end(),
            thrust::make_permutation_iterator
            (
                k.getField().begin(),
                patch().faceCells().end()
            ),
            muw.end(),
            rhow.begin()
        )),
        mutw.begin(),
        mutkCalcMutFunctor(yPlusLam_,Cmu25,kappa_,E_)
    );

/*
    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus =
            Cmu25*y[faceI]*sqrt(k[faceCellI])/(muw[faceI]/rhow[faceI]);

        if (yPlus > yPlusLam_)
        {
            mutw[faceI] = muw[faceI]*(yPlus*kappa_/log(E_*yPlus) - 1);
        }
    }
*/
    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf
)
:
    mutWallFunctionFvPatchScalarField(wfpsf)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalargpuField> mutkWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalargpuField kwc(k.boundaryField()[patchi].patchInternalField());
    const scalargpuField& muw = turbModel.mu().boundaryField()[patchi];
    const scalargpuField& rhow = turbModel.rho().boundaryField()[patchi];

    return pow025(Cmu_)*y*sqrt(kwc)/(muw/rhow);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutkWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
