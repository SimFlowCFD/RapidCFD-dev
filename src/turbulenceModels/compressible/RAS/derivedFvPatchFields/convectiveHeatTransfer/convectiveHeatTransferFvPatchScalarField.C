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

#include "convectiveHeatTransferFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    L_(1.0)
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    L_(ptf.L_)
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    L_(readScalar(dict.lookup("L")))
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& htcpsf
)
:
    fixedValueFvPatchScalarField(htcpsf),
    L_(htcpsf.L_)
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& htcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(htcpsf, iF),
    L_(htcpsf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

struct htcFunctor{
	const scalar L;
	htcFunctor(const scalar L_):L(L_){}
	__HOST____DEVICE__
	scalar operator () (const scalar& rhow, const thrust::tuple<vector,vector,scalar,scalar,scalar>& t){
		vector Uc = thrust::get<0>(t);
		vector Uw = thrust::get<1>(t);
		scalar muw = thrust::get<2>(t);
		scalar Pr = thrust::get<3>(t);
		scalar kappaw = thrust::get<4>(t);
		
		scalar Re = rhow*mag(Uc - Uw)*L/muw;

        if (Re < 5.0E+05)
        {
            return 0.664*sqrt(Re)*cbrt(Pr)*kappaw/L;
        }
        else
        {
            return 0.037*pow(Re, 0.8)*cbrt(Pr)*kappaw/L;
        }
	}
};

void convectiveHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalargpuField alphaEffw(turbModel.alphaEff(patchi));
    const scalargpuField& muw = turbModel.mu().boundaryField()[patchi];
    const scalargpuField& rhow = turbModel.rho().boundaryField()[patchi];
    const vectorgpuField& Uc = turbModel.U();
    const vectorgpuField& Uw = turbModel.U().boundaryField()[patchi];
    const scalargpuField& Tw = turbModel.thermo().T().boundaryField()[patchi];
    const scalargpuField& pw = turbModel.thermo().p().boundaryField()[patchi];
    const scalargpuField Cpw(turbModel.thermo().Cp(pw, Tw, patchi));

    const scalargpuField kappaw(Cpw*alphaEffw);
    const scalargpuField Pr(muw*Cpw/kappaw);

    scalargpuField& htc = *this;
    
    thrust::transform(rhow.begin(),rhow.end(),
                      thrust::make_zip_iterator(thrust::make_tuple(
                                                thrust::make_permutation_iterator(Uc.begin(),patch().faceCells().begin()),
                                                Uw.begin(),
                                                muw.begin(),
                                                Pr.begin(),
                                                kappaw.begin()
                                                )),
                      htc.begin(),
                      htcFunctor(L_));
    
/*
    forAll(htc, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar Re = rhow[faceI]*mag(Uc[faceCellI] - Uw[faceI])*L_/muw[faceI];

        if (Re < 5.0E+05)
        {
            htc[faceI] = 0.664*sqrt(Re)*cbrt(Pr[faceI])*kappaw[faceI]/L_;
        }
        else
        {
            htc[faceI] = 0.037*pow(Re, 0.8)*cbrt(Pr[faceI])*kappaw[faceI]/L_;
        }
    }
*/
    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void convectiveHeatTransferFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    convectiveHeatTransferFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
