/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "epsilonLowReWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline scalar yPlusLam
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

struct EpsilonLowReCalculateEpsilonFunctor : public std::unary_function<label,scalar>
{
    const scalar yPlusLam;
    const scalar Cmu25;
    const scalar Cmu75;
    const scalar kappa;
    const scalar* cornerWeights;
    const scalar* y;
    const scalar* k;
    const scalar* nuw;

    EpsilonLowReCalculateEpsilonFunctor
    (
        const scalar yPlusLam_,
        const scalar Cmu25_,
        const scalar Cmu75_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* nuw_
    ):
        yPlusLam(yPlusLam_),
        Cmu25(Cmu25_),
        Cmu75(Cmu75_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        nuw(nuw_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI,const label& faceI)
    {
        scalar yPlus = Cmu25*sqrt(k[cellI])*y[faceI]/nuw[faceI];

        scalar w = cornerWeights[faceI];

        if (yPlus > yPlusLam)
        {
            return w*Cmu75*pow(k[cellI], 1.5)/(kappa*y[faceI]);
        }
        else
        {
            return w*2.0*k[cellI]*nuw[faceI]/sqr(y[faceI]);
        }
    }
};

struct EpsilonLowReCalculateGFunctor : public std::unary_function<label,scalar>
{
	const scalar Cmu25;
	const scalar kappa;
	const scalar* cornerWeights;
	const scalar* y;
	const scalar* k;
	const scalar* nuw;
	const scalar* nutw;
	const scalar* magGradUw;

    EpsilonLowReCalculateGFunctor
    (
        const scalar Cmu25_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* nuw_,
        const scalar* nutw_,
        const scalar* magGradUw_
    ):
        Cmu25(Cmu25_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        nuw(nuw_),
        nutw(nutw_),
        magGradUw(magGradUw_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        return
            w
           *(nutw[faceI] + nuw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[cellI])
           /(kappa*y[faceI]);
    }
};

void epsilonLowReWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const gpuList<scalar>& cornerWeights,
    const fvPatch& patch,
    scalargpuField& G,
    scalargpuField& epsilon
)
{
    const label patchI = patch.index();

    const scalargpuField& y = turbulence.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tnu = turbulence.nu();
    const scalargpuField& nuw = tnu().boundaryField()[patchI];

    const tmp<volScalarField> tnut = turbulence.nut();
    const volScalarField& nut = tnut();
    const scalargpuField& nutw = nut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalargpuField magGradUw(mag(Uw.snGrad()));

    // Set epsilon and G
                                               
    matrixPatchOperation
    (
        patchI,
        epsilon,
        patch.boundaryMesh().mesh().lduAddr(),
        EpsilonLowReCalculateEpsilonFunctor
        (
            yPlusLam_,
            Cmu25,
            Cmu75,
            kappa_,
            cornerWeights.data(),
            y.data(),
            k.getField().data(),
            nuw.data()
        )
    );
	
    matrixPatchOperation
    (
        patchI,
        G,
        patch.boundaryMesh().mesh().lduAddr(),
        EpsilonLowReCalculateGFunctor
        (
            Cmu25,
            kappa_,
            cornerWeights.data(),
            y.data(),
            k.getField().data(),
            nuw.data(),
            nutw.data(),
            magGradUw.data()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF),
    yPlusLam_(yPlusLam(kappa_, E_))
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    yPlusLam_(ptf.yPlusLam_)
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF, dict),
    yPlusLam_(yPlusLam(kappa_, E_))
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf),
    yPlusLam_(ewfpsf.yPlusLam_)
{}


epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf, iF),
    yPlusLam_(ewfpsf.yPlusLam_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
