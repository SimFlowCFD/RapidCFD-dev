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

#include "epsilonLowReWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
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
    const scalar* muw;
    const scalar* rhow;

    EpsilonLowReCalculateEpsilonFunctor
    (
        const scalar yPlusLam_,
        const scalar Cmu25_,
        const scalar Cmu75_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* muw_,
        const scalar* rhow_
    ):
        yPlusLam(yPlusLam_),
        Cmu25(Cmu25_),
        Cmu75(Cmu75_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        muw(muw_),
        rhow(rhow_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar yPlus = Cmu25*sqrt(k[cellI])*y[faceI]/(muw[faceI]/rhow[faceI]);

        scalar w = cornerWeights[faceI];

        if (yPlus > yPlusLam)
        {
            return w*Cmu75*pow(k[cellI], 1.5)/(kappa*y[faceI]);
        }
        else
        {
            return w*2.0*k[cellI]*muw[faceI]/rhow[faceI]/sqr(y[faceI]);
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
    const scalar* muw;
    const scalar* mutw;
    const scalar* magGradUw;

    EpsilonLowReCalculateGFunctor
    (
        const scalar Cmu25_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* muw_,
        const scalar* mutw_,
        const scalar* magGradUw_
    ):
        Cmu25(Cmu25_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        muw(muw_),
        mutw(mutw_),
        magGradUw(magGradUw_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        return
            w
           *(mutw[faceI] + muw[faceI])
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

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalargpuField& muw = tmu().boundaryField()[patchI];

    const tmp<volScalarField> tmut = turbulence.mut();
    const volScalarField& mut = tmut();
    const scalargpuField& mutw = mut.boundaryField()[patchI];

    const scalargpuField& rhow = turbulence.rho().boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalargpuField magGradUw(mag(Uw.snGrad()));
    
	
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
            muw.data(),
            rhow.data()
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
            muw.data(),
            mutw.data(),
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

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
