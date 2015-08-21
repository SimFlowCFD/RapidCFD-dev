/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
namespace Foam
{
    template<class Limiter,class phiType, class gradPhiType>
    struct LimitedSchemeCalcLimiterFunctor
    {
        const Limiter limiter;

        LimitedSchemeCalcLimiterFunctor(const Limiter& _limiter):limiter(_limiter){}

        template<class Tuple>
        __HOST____DEVICE__
        scalar operator()(const scalar& w, const Tuple& t)
	{
            return limiter.limiter
            (
                w,
                thrust::get<0>(t),
                thrust::get<1>(t),
                thrust::get<2>(t),
                thrust::get<3>(t),
                thrust::get<4>(t),
                thrust::get<5>(t) - thrust::get<6>(t)
            );
        }
    };
}

template<class Type, class Limiter, template<class> class LimitFunc>
void Foam::LimitedScheme<Type, Limiter, LimitFunc>::calcLimiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    surfaceScalarField& limiterField
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh> >
        tlPhi = LimitFunc<Type>()(phi);

    const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
        lPhi = tlPhi();

    tmp<GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh> >
        tgradc(fvc::grad(lPhi));
    const GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>&
        gradc = tgradc();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const labelgpuList& owner = mesh.owner();
    const labelgpuList& neighbour = mesh.neighbour();

    const vectorgpuField& C = mesh.C();

    scalargpuField& pLim = limiterField.internalField();
    
    thrust::transform
    (
        CDweights.getField().begin(),
        CDweights.getField().end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            this->faceFlux_.getField().begin(),
            thrust::make_permutation_iterator
            (
                lPhi.getField().begin(),
                owner.begin()
            ),
            thrust::make_permutation_iterator
            (
                lPhi.getField().begin(),
                neighbour.begin()
            ),
            thrust::make_permutation_iterator
            (
                gradc.getField().begin(),
                owner.begin()
            ),
            thrust::make_permutation_iterator
            (
                gradc.getField().begin(),
                neighbour.begin()
            ),
            thrust::make_permutation_iterator
            (
                C.begin(),
                neighbour.begin()
            ),
            thrust::make_permutation_iterator
            (
                C.begin(),
                owner.begin()
            )
        )),
        pLim.begin(),
        LimitedSchemeCalcLimiterFunctor
        <
            Limiter,
            typename Limiter::phiType,
            typename Limiter::gradPhiType
        >
        (
            static_cast<const Limiter&>(*this)
        )
    );

    surfaceScalarField::GeometricBoundaryField& bLim =
        limiterField.boundaryField();

    forAll(bLim, patchi)
    {
        scalargpuField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalargpuField& pCDweights = CDweights.boundaryField()[patchi];
            const scalargpuField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            const gpuField<typename Limiter::phiType> plPhiP
            (
                lPhi.boundaryField()[patchi].patchInternalField()
            );
            const gpuField<typename Limiter::phiType> plPhiN
            (
                lPhi.boundaryField()[patchi].patchNeighbourField()
            );
            const gpuField<typename Limiter::gradPhiType> pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );
            const gpuField<typename Limiter::gradPhiType> pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorgpuField pd(CDweights.boundaryField()[patchi].patch().delta());

            thrust::transform
            (
                pCDweights.begin(),
                pCDweights.end(),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    pFaceFlux.begin(),
                    plPhiP.begin(),
                    plPhiN.begin(),
                    pGradcP.begin(),
                    pGradcN.begin(),
                    pd.begin(),
                    thrust::make_constant_iterator(vector(0,0,0))
                )),
                pLim.begin(),
                LimitedSchemeCalcLimiterFunctor
                <
                    Limiter,
                    typename Limiter::phiType,
                    typename Limiter::gradPhiType
                >
                (
                    static_cast<const Limiter&>(*this)
                )
            );
        }
        else
        {
            pLim = 1.0;
        }
    }
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
Foam::tmp<Foam::surfaceScalarField>
Foam::LimitedScheme<Type, Limiter, LimitFunc>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    const word limiterFieldName(type() + "Limiter(" + phi.name() + ')');

    if (this->mesh().cache("limiter"))
    {
        if (!mesh.foundObject<surfaceScalarField>(limiterFieldName))
        {
            surfaceScalarField* limiterField
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        limiterFieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimless
                )
            );
            
            const polyMesh& pm = mesh;
            pm.store(limiterField);
        }

        surfaceScalarField& limiterField =
            const_cast<surfaceScalarField&>
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    limiterFieldName
                )
            );

        calcLimiter(phi, limiterField);

        return limiterField;
    }
    else
    {
        tmp<surfaceScalarField> tlimiterField
        (
            new surfaceScalarField
            (
                IOobject
                (
                    limiterFieldName,
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimless
            )
        );

        calcLimiter(phi, tlimiterField());

        return tlimiterField;
    }
}


// ************************************************************************* //
