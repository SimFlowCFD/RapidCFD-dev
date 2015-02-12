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

#include "linearUpwindV.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class ProdType>
struct linearUpwindVSfCorrFunctor
{
    const Type zero;

    linearUpwindVSfCorrFunctor
    (
    ):
       zero(pTraits<Type>::zero)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    Type operator()(const scalar& faceFlux, const Tuple& t)
    {
        const scalar& w = thrust::get<0>(t);
        const Type& vfn = thrust::get<1>(t);
        const Type& vfo = thrust::get<2>(t);
        const vector& Cf = thrust::get<3>(t);
        const vector& Co = thrust::get<4>(t);
        const vector& Cn = thrust::get<5>(t);
        const ProdType& gradVfo = thrust::get<6>(t);
        const ProdType& gradVfn = thrust::get<7>(t);

        Type maxCorr;
        Type sfCorr;

        if (faceFlux > 0.0)
        {
            maxCorr = (1.0 - w)*(vfn - vfo);

            sfCorr = (Cf - Co) & gradVfo;
        }
        else
        {
            maxCorr = w*(vfo - vfn);

            sfCorr = (Cf - Cn) & gradVfn;
        }

        scalar sfCorrs = magSqr(sfCorr);
        scalar maxCorrs = sfCorr & maxCorr;

        if (sfCorrs > 0)
        {
            if (maxCorrs < 0)
            {
                return zero;
            }
            else if (sfCorrs > maxCorrs)
            {
                return sfCorr * maxCorrs/(sfCorrs + VSMALL);
            }
        }
        else if (sfCorrs < 0)
        {
            if (maxCorrs > 0)
            {
                return zero;
            }
            else if (sfCorrs < maxCorrs)
            {
                return sfCorr * maxCorrs/(sfCorrs - VSMALL);
            }
        }
    }
};


template<class Type, class ProdType>
struct linearUpwindVPatchSfCorrFunctor
{
    const Type zero;

    linearUpwindVPatchSfCorrFunctor
    (
    ):
       zero(pTraits<Type>::zero)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    Type operator()(const scalar& pFaceFlux, const Tuple& t)
    {
        const scalar& pW = thrust::get<0>(t);
        const Type& vfo = thrust::get<1>(t);
        const Type& pVfNei = thrust::get<2>(t);
        const vector& pCf = thrust::get<3>(t);
        const vector& Co = thrust::get<4>(t);
        const vector& pd = thrust::get<5>(t);
        const ProdType& gradVfo = thrust::get<6>(t);
        const ProdType& pGradVfNei = thrust::get<7>(t);

        vector maxCorr;
        Type pSfCorr;

        if (pFaceFlux > 0)
        {
            pSfCorr = (pCf - Co) & gradVfo;

            maxCorr = (1.0 - pW)*(pVfNei - vfo);
        }
        else
        {
            pSfCorr = (pCf - pd - Co) & pGradVfNei;

            maxCorr = pW*(vfo - pVfNei);
        }

        scalar pSfCorrs = magSqr(pSfCorr);
        scalar maxCorrs = pSfCorr & maxCorr;

        if (pSfCorrs > 0)
        {
            if (maxCorrs < 0)
            {
                return zero;
            }
            else if (pSfCorrs > maxCorrs)
            {
                return pSfCorr * maxCorrs/(pSfCorrs + VSMALL);
            }
        }
        else if (pSfCorrs < 0)
        {
            if (maxCorrs > 0)
            {
                return zero;
            }
            else if (pSfCorrs < maxCorrs)
            {
                return pSfCorr * maxCorrs/(pSfCorrs - VSMALL);
            }
        }
    }
};

}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::linearUpwindV<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "linearUpwindV::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>
            (
                vf.name(),
                vf.dimensions(),
                pTraits<Type>::zero
            )
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

    const surfaceScalarField& faceFlux = this->faceFlux_;
    const surfaceScalarField& w = mesh.weights();

    const labelgpuList& own = mesh.owner();
    const labelgpuList& nei = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > tgradVf = gradScheme_().grad(vf, gradSchemeName_);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gradVf = tgradVf();

    thrust::transform
    (
        faceFlux.getField().begin(),
        faceFlux.getField().end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            w.getField().begin(), 
            thrust::make_permutation_iterator
            (
                vf.getField().begin(),
                nei.begin()
            ),
            thrust::make_permutation_iterator
            (
                vf.getField().begin(),
                own.begin()
            ),
            Cf.getField().begin(),
            thrust::make_permutation_iterator
            (
                C.getField().begin(),
                own.begin()
            ),
            thrust::make_permutation_iterator
            (
                C.getField().begin(),
                nei.begin()
            ),
            thrust::make_permutation_iterator
            (
                gradVf.getField().begin(),
                own.begin()
            ),
            thrust::make_permutation_iterator
            (
                gradVf.getField().begin(),
                nei.begin()
            )
        )),
        sfCorr.getField().begin(),
        linearUpwindVSfCorrFunctor
        <
            Type,
            typename outerProduct<vector, Type>::type
        >()
    );

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        GeometricBoundaryField& bSfCorr = sfCorr.boundaryField();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelgpuList& pOwner =
                mesh.boundary()[patchi].faceCells();

            const vectorgpuField& pCf = Cf.boundaryField()[patchi];
            const scalargpuField& pW = w.boundaryField()[patchi];

            const scalargpuField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const gpuField<typename outerProduct<vector, Type>::type> pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            const gpuField<Type> pVfNei
            (
                vf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorgpuField pd(Cf.boundaryField()[patchi].patch().delta());

            thrust::transform
            (
                pFaceFlux.begin(),
                pFaceFlux.end(),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    pW.begin(),
                    thrust::make_permutation_iterator
                    (
                        vf.getField().begin(),
                        pOwner.begin()
                    ),
                    pVfNei.begin(),
                    pCf.begin(),
                    thrust::make_permutation_iterator
                    (
                        C.getField().begin(),
                        pOwner.begin()
                    ),
                    pd.begin(),
                    thrust::make_permutation_iterator
                    (
                        gradVf.getField().begin(),
                        pOwner.begin()
                    ),
                    pGradVfNei.begin()
                )),
                pSfCorr.begin(), 
                linearUpwindVPatchSfCorrFunctor
                <
                    Type,
                    typename outerProduct<vector, Type>::type
                >()
            );
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makelimitedSurfaceInterpolationTypeScheme(linearUpwindV, vector)
}

// ************************************************************************* //
