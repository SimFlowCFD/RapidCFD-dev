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

#include "linearUpwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type,class ProdType>
struct linearUpwindCorrectionFunctor: std::binary_function<label,thrust::tuple<scalar,vector>,Type>{
	const label* own;
	const label* nei;
	const vector* C;
	const ProdType* gradf;
	linearUpwindCorrectionFunctor(const label* _own,const label* _nei, const vector* _C, const ProdType* _gradf):
	                              own(_own),nei(_nei),C(_C),gradf(_gradf){}
    __HOST____DEVICE__
    Type operator()(const label& id,const thrust::tuple<scalar,vector>& t){
        label celli = thrust::get<0>(t) > 0?own[id]:nei[id];
        return (thrust::get<1>(t) - C[celli]) & gradf[celli];
    }
};

template<class Type,class ProdType>
struct linearUpwindPatchCorrectionFunctor{
    __HOST____DEVICE__
    Type operator()(const thrust::tuple<scalar,vector,vector,vector,ProdType,ProdType>& t){
        if(thrust::get<0>(t) > 0)
        {
             return (thrust::get<1>(t) - thrust::get<2>(t)) & thrust::get<4>(t);
        }
        else
        {
             return (thrust::get<1>(t) - thrust::get<2>(t) - thrust::get<3>(t)) & thrust::get<5>(t);
        }
    }
};
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::linearUpwind<Type>::correction
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
                "linearUpwind::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelgpuList& owner = mesh.owner();
    const labelgpuList& neighbour = mesh.neighbour();

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
/*
    forAll(faceFlux, facei)
    {
        label celli = (faceFlux[facei] > 0) ? owner[facei] : neighbour[facei];
        sfCorr[facei] = (Cf[facei] - C[celli]) & gradVf[celli];
    }
*/
    thrust::transform(thrust::make_counting_iterator(0),thrust::make_counting_iterator(0)+faceFlux.size(),
                      thrust::make_zip_iterator(thrust::make_tuple(faceFlux.getField().begin(),
                                                                   Cf.getField().begin()
                                                                  )
                                               ),
                     sfCorr.getField().begin(),
                     linearUpwindCorrectionFunctor<Type,typename outerProduct<vector, Type>::type>(
                                                   owner.data(),
                                                   neighbour.data(),
                                                   C.getField().data(),
                                                   gradVf.getField().data()
                                                   ));


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

            const scalargpuField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const gpuField<typename outerProduct<vector, Type>::type> pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorgpuField pd(Cf.boundaryField()[patchi].patch().delta());
/*
            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
                }
                else
                {
                    pSfCorr[facei] =
                        (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];
                }
            }
*/
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(pFaceFlux.begin(),
                                                                   pCf.begin(),
                                                                   thrust::make_permutation_iterator(C.getField().begin(),pOwner.begin()),
                                                                   pd.begin(),
                                                                   thrust::make_permutation_iterator(gradVf.getField().begin(),pOwner.begin()),
                                                                   pGradVfNei.begin()
                                                                  )
                                               ),
                      thrust::make_zip_iterator(thrust::make_tuple(pFaceFlux.end(),
                                                                   pCf.end(),
                                                                   thrust::make_permutation_iterator(C.getField().begin(),pOwner.end()),
                                                                   pd.end(),
                                                                   thrust::make_permutation_iterator(gradVf.getField().end(),pOwner.begin()),
                                                                   pGradVfNei.end()
                                                                  )
                                               ),
                     pSfCorr.begin(),
                     linearUpwindPatchCorrectionFunctor<Type,typename outerProduct<vector, Type>::type>());

        }
    }

    return tsfCorr;
}


namespace Foam
{
    //makelimitedSurfaceInterpolationScheme(linearUpwind)
    makelimitedSurfaceInterpolationTypeScheme(linearUpwind, scalar)
    makelimitedSurfaceInterpolationTypeScheme(linearUpwind, vector)
}

// ************************************************************************* //
