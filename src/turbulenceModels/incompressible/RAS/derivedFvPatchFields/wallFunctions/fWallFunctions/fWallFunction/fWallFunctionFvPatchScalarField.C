/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "fWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "v2f.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void fWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("fWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void fWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


scalar fWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fWallFunctionFvPatchScalarField& v2wfpsf
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf),
    Cmu_(v2wfpsf.Cmu_),
    kappa_(v2wfpsf.kappa_),
    E_(v2wfpsf.E_),
    yPlusLam_(v2wfpsf.yPlusLam_)
{
    checkType();
}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fWallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf, iF),
    Cmu_(v2wfpsf.Cmu_),
    kappa_(v2wfpsf.kappa_),
    E_(v2wfpsf.E_),
    yPlusLam_(v2wfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

struct fWallFunctionFunctor
{
    const scalar Cmu25_;
    const scalar yPlusLam_;
 
    fWallFunctionFunctor(scalar Cmu25,scalar yPlusLam):
                         Cmu25_(Cmu25),yPlusLam_(yPlusLam){}

    __HOST____DEVICE__
    scalar operator()(const scalar& y, const thrust::tuple<scalar,scalar,scalar,scalar>& t)
    {
        const scalar& kc = thrust::get<0>(t);
        const scalar& nuw = thrust::get<1>(t);
        const scalar& v2c = thrust::get<2>(t);
        const scalar& epsc = thrust::get<3>(t);

        scalar uTau = Cmu25_*sqrt(kc);

        scalar yPlus = uTau*y/nuw;

        if (yPlus > yPlusLam_)
        {
            scalar N = 6.0;

            return (N*v2c*epsc/(sqr(kc) + ROOTVSMALL))
                   / (sqr(uTau) + ROOTVSMALL);
        }
        else
        {
            return 0.0;
        }
    }
};

void fWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const v2f& v2fModel = refCast<const v2f>(turbulence);

    const scalargpuField& y = v2fModel.y()[patchI];

    const tmp<volScalarField> tk = v2fModel.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tepsilon = v2fModel.epsilon();
    const volScalarField& epsilon = tepsilon();

    const tmp<volScalarField> tv2 = v2fModel.v2();
    const volScalarField& v2 = tv2();

    const tmp<volScalarField> tnu = v2fModel.nu();
    const scalargpuField& nuw = tnu().boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    scalargpuField& f = *this;

    thrust::transform
    (
        y.begin(),
        y.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator
            (
                k.getField().begin(),
                patch().faceCells().begin()
            ),
            nuw.begin(),
            thrust::make_permutation_iterator
            (
                v2.getField().begin(),
                patch().faceCells().begin()
            ),
            thrust::make_permutation_iterator
            (
                epsilon.getField().begin(),
                patch().faceCells().begin()
            )
        )),
        f.begin(),
        fWallFunctionFunctor(Cmu25,yPlusLam_)
    );


    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void fWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void fWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
