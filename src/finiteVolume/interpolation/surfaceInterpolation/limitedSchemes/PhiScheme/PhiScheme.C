/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	template<class Type,class PhiLimiter>
	struct PhiSchemeLimiterFunctor{
		const PhiLimiter& limiter;
		PhiSchemeLimiterFunctor(const PhiLimiter& _limiter): limiter(_limiter) {}
		__HOST____DEVICE__
		scalar operator ()(const thrust::tuple<scalar,scalar,Type,Type,vector,scalar>& t){
			return limiter.limiter(thrust::get<0>(t),
			                       thrust::get<1>(t),
			                       thrust::get<2>(t),
			                       thrust::get<3>(t),
			                       thrust::get<4>(t),
			                       thrust::get<5>(t));
		}
	};
}

template<class Type, class PhiLimiter>
Foam::tmp<Foam::surfaceScalarField>
Foam::PhiScheme<Type, PhiLimiter>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                "PhiLimiter",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );
    surfaceScalarField& Limiter = tLimiter();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const labelgpuList& owner = mesh.owner();
    const labelgpuList& neighbour = mesh.neighbour();

    tmp<surfaceScalarField> tUflux = this->faceFlux_;

    if (this->faceFlux_.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const volScalarField& rho =
            phi.db().objectRegistry::template lookupObject<volScalarField>
            ("rho");

        tUflux = this->faceFlux_/fvc::interpolate(rho);
    }
    else if (this->faceFlux_.dimensions() != dimVelocity*dimArea)
    {
        FatalErrorIn
        (
            "PhiScheme<PhiLimiter>::limiter"
            "(const GeometricField<Type, fvPatchField, volMesh>& phi)"
        )   << "dimensions of faceFlux are not correct"
            << exit(FatalError);
    }

    const surfaceScalarField& Uflux = tUflux();

    scalargpuField& pLimiter = Limiter.internalField();
/*
    forAll(pLimiter, face)
    {
        pLimiter[face] = PhiLimiter::limiter
        (
            CDweights[face],
            Uflux[face],
            phi[owner[face]],
            phi[neighbour[face]],
            Sf[face],
            magSf[face]
        );
    }
    */
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(CDweights.getField().begin(),
                                                                   Uflux.getField().begin(),
                                                                   thrust::make_permutation_iterator(phi.getField().begin(),owner.begin()),
                                                                   thrust::make_permutation_iterator(phi.getField().begin(),neighbour.begin()),
                                                                   Sf.getField().begin(),
                                                                   magSf.getField().begin())),
                     thrust::make_zip_iterator(thrust::make_tuple(CDweights.getField().end(),
                                                                   Uflux.getField().end(),
                                                                   thrust::make_permutation_iterator(phi.getField().begin(),owner.end()),
                                                                   thrust::make_permutation_iterator(phi.getField().begin(),neighbour.end()),
                                                                   Sf.getField().end(),
                                                                   magSf.getField().end())),
                      pLimiter.begin(),
                      PhiSchemeLimiterFunctor<Type,PhiLimiter>(*this));


    surfaceScalarField::GeometricBoundaryField& bLimiter =
        Limiter.boundaryField();

    forAll(bLimiter, patchI)
    {
        scalargpuField& pLimiter = bLimiter[patchI];

        if (bLimiter[patchI].coupled())
        {
            const scalargpuField& pCDweights = CDweights.boundaryField()[patchI];
            const vectorgpuField& pSf = Sf.boundaryField()[patchI];
            const scalargpuField& pmagSf = magSf.boundaryField()[patchI];
            const scalargpuField& pFaceFlux = Uflux.boundaryField()[patchI];

            const gpuField<Type> pphiP
            (
                phi.boundaryField()[patchI].patchInternalField()
            );
            const gpuField<Type> pphiN
            (
                phi.boundaryField()[patchI].patchNeighbourField()
            );
/*
            forAll(pLimiter, face)
            {
                pLimiter[face] = PhiLimiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    pphiP[face],
                    pphiN[face],
                    pSf[face],
                    pmagSf[face]
                );
            }
            * 
            */
            thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(pCDweights.begin(),
                                                                   pFaceFlux.begin(),
                                                                   pphiP.begin(),
                                                                   pphiN.begin(),
                                                                   pSf.begin(),
                                                                   pmagSf.begin())),
                               thrust::make_zip_iterator(thrust::make_tuple(pCDweights.end(),
                                                                   pFaceFlux.end(),
                                                                   pphiP.end(),
                                                                   pphiN.end(),
                                                                   pSf.end(),
                                                                   pmagSf.end())),
                                pLimiter.begin(),
                                   PhiSchemeLimiterFunctor<Type,PhiLimiter>(*this));
        }
        else
        {
            pLimiter = 1.0;
        }
    }

    return tLimiter;
}


// ************************************************************************* //
