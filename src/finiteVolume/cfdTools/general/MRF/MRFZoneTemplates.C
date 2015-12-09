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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    template<bool relative>
    struct MFRZoneRhoFluxFunctor
    {
        const vector omega;
        const vector origin;

        MFRZoneRhoFluxFunctor
        (
            vector _omega, 
            vector _origin
        ): 
            omega(_omega), 
            origin(_origin)
        {}

        __host__ __device__                                        //rho   Cfi     Sfi
        scalar operator () (const scalar& phi, const thrust::tuple<scalar,vector,vector>& t) const
        {
            scalar delta = thrust::get<0>(t)*(omega ^ (thrust::get<1>(t) - origin)) & thrust::get<2>(t);
            if(relative)
                return phi - delta;
            else
                return phi + delta;
        }
    };
}

template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = omega_->value(mesh_.time().timeOutputValue())*axis_;

    const vectorgpuField& Cfi = Cf.internalField();
    const vectorgpuField& Sfi = Sf.internalField();
    scalargpuField& phii = phi.internalField();

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            phii.begin(),
            internalFaces_.begin()
        ),
        thrust::make_permutation_iterator
        (
            phii.begin(),
            internalFaces_.end()
        ),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator(rho.getField().begin(), internalFaces_.begin()),
            thrust::make_permutation_iterator(Cfi.begin(), internalFaces_.begin()),
            thrust::make_permutation_iterator(Sfi.begin(), internalFaces_.begin())
        )),
        thrust::make_permutation_iterator
        (
            phii.begin(),
            internalFaces_.begin()
        ),
        MFRZoneRhoFluxFunctor<true>(Omega,origin_)
    );

    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryField());
}


template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = omega_->value(mesh_.time().timeOutputValue())*axis_;

    // Included patches
    forAll(includedFaces_, patchi)
    {
        const labelgpuList& faces = includedFaces_[patchi];
        scalargpuField& phii = phi[patchi];
		
        thrust::fill
        (
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            thrust::make_permutation_iterator(phii.begin(),faces.end()),
            scalar(0.0)
        );
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        const labelgpuList& faces = excludedFaces_[patchi];
        scalargpuField& phii = phi[patchi];
		
        thrust::transform
        (
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            thrust::make_permutation_iterator(phii.begin(),faces.end()),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator(rho[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Cf.boundaryField()[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Sf.boundaryField()[patchi].begin(),faces.begin())
            )),
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            MFRZoneRhoFluxFunctor<true>(Omega,origin_)
        );
    }
}


template<class RhoFieldType>
void Foam::MRFZone::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = omega_->value(mesh_.time().timeOutputValue())*axis_;

    const vectorgpuField& Cfi = Cf.internalField();
    const vectorgpuField& Sfi = Sf.internalField();
    scalargpuField& phii = phi.internalField();  
    
    thrust::transform
    (
        thrust::make_permutation_iterator(phii.begin(),internalFaces_.begin()),
        thrust::make_permutation_iterator(phii.begin(),internalFaces_.end()),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator(rho.getField().begin(),internalFaces_.begin()),
            thrust::make_permutation_iterator(Cfi.begin(),internalFaces_.begin()),
            thrust::make_permutation_iterator(Sfi.begin(),internalFaces_.begin())
        )),
        thrust::make_permutation_iterator(phii.begin(),internalFaces_.begin()),
        MFRZoneRhoFluxFunctor<false>(Omega,origin_)
    );

    // Included patches
    forAll(includedFaces_, patchi)
    {
        const labelgpuList& faces = includedFaces_[patchi];
        scalargpuField& phii = phi.boundaryField()[patchi];
		
        thrust::transform
        (
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            thrust::make_permutation_iterator(phii.begin(),faces.end()),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator(rho.boundaryField()[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Cf.boundaryField()[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Sf.boundaryField()[patchi].begin(),faces.begin())
            )),
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            MFRZoneRhoFluxFunctor<false>(Omega,origin_)
        );
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        const labelgpuList& faces = excludedFaces_[patchi];
        scalargpuField& phii = phi.boundaryField()[patchi];
		
        thrust::transform
        (
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            thrust::make_permutation_iterator(phii.begin(),faces.end()),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator(rho.boundaryField()[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Cf.boundaryField()[patchi].begin(),faces.begin()),
                thrust::make_permutation_iterator(Sf.boundaryField()[patchi].begin(),faces.begin())
            )),
            thrust::make_permutation_iterator(phii.begin(),faces.begin()),
            MFRZoneRhoFluxFunctor<false>(Omega,origin_)
        );
    }
}


// ************************************************************************* //
