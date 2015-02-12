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

#include "leastSquaresGrad.H"
#include "leastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class Type,class GradType,bool add>
struct leastSquaresCalcGradFunctor{
    __HOST____DEVICE__
    GradType operator()(const GradType& g,const thrust::tuple<vector,Type,Type>& t){
        GradType d = thrust::get<0>(t)*(thrust::get<1>(t) - thrust::get<2>(t));
        if(add)
        {
             return g + d;
        } 
        else
        {
             return g - d;
        }
    }
};
}

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::leastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();

    // Get reference to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelgpuList& own = mesh.owner();
    const labelgpuList& nei = mesh.neighbour();
/*
    forAll(own, facei)
    {
        register label ownFaceI = own[facei];
        register label neiFaceI = nei[facei];

        Type deltaVsf = vsf[neiFaceI] - vsf[ownFaceI];

        lsGrad[ownFaceI] += ownLs[facei]*deltaVsf;
        lsGrad[neiFaceI] -= neiLs[facei]*deltaVsf;
    }
*/
    thrust::transform(thrust::make_permutation_iterator(lsGrad.getField().begin(),own.begin()),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),own.end()),
                      thrust::make_zip_iterator(thrust::make_tuple(ownLs.getField().begin(),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),nei.begin()),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),own.begin())
                                                                  )),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),own.begin()),
                      leastSquaresCalcGradFunctor<Type,GradType,true>());

    thrust::transform(thrust::make_permutation_iterator(lsGrad.getField().begin(),nei.begin()),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),nei.end()),
                      thrust::make_zip_iterator(thrust::make_tuple(ownLs.getField().begin(),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),nei.begin()),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),own.begin())
                                                                  )),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),nei.begin()),
                      leastSquaresCalcGradFunctor<Type,GradType,false>());

    // Boundary faces
    forAll(vsf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const labelgpuList& faceCells =
            lsGrad.boundaryField()[patchi].patch().faceCells();

        if (vsf.boundaryField()[patchi].coupled())
        {
            const gpuField<Type> neiVsf
            (
                vsf.boundaryField()[patchi].patchNeighbourField()
            );
/*
            forAll(neiVsf, patchFaceI)
            {
                lsGrad[faceCells[patchFaceI]] +=
                    patchOwnLs[patchFaceI]
                   *(neiVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
            }
*/
    thrust::transform(thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.begin()),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.end()),
                      thrust::make_zip_iterator(thrust::make_tuple(patchOwnLs.getField().begin(),
                                                                   neiVsf.begin(),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),faceCells.begin())
                                                                  )),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.begin()),
                      leastSquaresCalcGradFunctor<Type,GradType,true>());
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
/*
            forAll(patchVsf, patchFaceI)
            {
                lsGrad[faceCells[patchFaceI]] +=
                     patchOwnLs[patchFaceI]
                    *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
            }
*/
    thrust::transform(thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.begin()),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.end()),
                      thrust::make_zip_iterator(thrust::make_tuple(patchOwnLs.getField().begin(),
                                                                   patchVsf.begin(),
                                                                   thrust::make_permutation_iterator(vsf.getField().begin(),faceCells.begin())
                                                                  )),
                      thrust::make_permutation_iterator(lsGrad.getField().begin(),faceCells.begin()),
                      leastSquaresCalcGradFunctor<Type,GradType,true>());
        }
    }


    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// ************************************************************************* //
