/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

struct reconstructMagFunctor
{
    const vector* Sf;
    const vector* Cf;
    const vector* C;
    const scalar* ssf;
    const scalar* magSf;

    reconstructMagFunctor
    (
        const vector* _Sf,
        const vector* _Cf,
        const vector* _C,
        const scalar* _ssf,
        const scalar* _magSf
    ):
        Sf(_Sf),
        Cf(_Cf),
        C(_C),
        ssf(_ssf),
        magSf(_magSf)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& celli, const label& facei)
    {
        return (Sf[facei] & (Cf[facei] - C[celli]))*ssf[facei]/magSf[facei];
    }   
};

struct reconstructMagPatchFunctor
{
    const vector* pSf;
    const vector* pCf;
    const vector* C;
    const scalar* psf;
    const scalar* pMagSf;

    reconstructMagPatchFunctor
    (
        const vector* _pSf,
        const vector* _pCf,
        const vector* _C,
        const scalar* _psf,
        const scalar* _pMagSf
    ):
        pSf(_pSf),
        pCf(_pCf),
        C(_C),
        psf(_psf),
        pMagSf(_pMagSf)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& celli, const label& pFacei)
    {
        return (pSf[pFacei] & (pCf[pFacei] - C[celli]))
               *psf[pFacei]/pMagSf[pFacei];
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> reconstructMag(const surfaceScalarField& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    tmp<volScalarField> treconField
    (
        new volScalarField
        (
            IOobject
            (
                "reconstruct("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "0",
                ssf.dimensions()/dimArea,
                scalar(0)
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalargpuField& rf = treconField();

    reconstructMagFunctor fun
    (
        Sf.getField().data(),
        Cf.getField().data(),
        C.getField().data(),
        ssf.getField().data(),
        magSf.getField().data()
    );

    matrixOperation
    (
        thrust::make_constant_iterator(scalar(0)),
        rf,
        mesh.lduAddr(),
        fun,
        fun,
        sumOp<scalar>(),
        minusOp<scalar>()
    );

    const surfaceScalarField::GeometricBoundaryField& bsf = ssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchScalarField& psf = bsf[patchi];

        const vectorgpuField& pCf = Cf.boundaryField()[patchi];
        const vectorgpuField& pSf = Sf.boundaryField()[patchi];
        const scalargpuField& pMagSf = magSf.boundaryField()[patchi];

        reconstructMagPatchFunctor pfun
        (
            pSf.data(),
            pCf.data(),
            C.getField().data(),
            psf.data(),
            pMagSf.data()
        );

        matrixPatchOperation
        (
            patchi,
            rf,
            mesh.lduAddr(),
            fun
        );
    }

    rf /= mesh.V();

    treconField().correctBoundaryConditions();

    return treconField;
}


tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>& tssf)
{
    tmp<volScalarField> tvf
    (
        fvc::reconstructMag(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
