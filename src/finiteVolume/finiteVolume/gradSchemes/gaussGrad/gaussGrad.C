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

#include "gaussGrad.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type,class GradType>
    struct gaussGradFunctor
    {
        const GradType zero;
        const vector* Sf;
        const Type* issf;
        const label* ownStart;
        const label* neiStart;
        const label* own;
        const label* nei;
        const label* losort;

        gaussGradFunctor
        (
            const GradType _zero,
            const vector* _Sf,
            const Type* _issf,
            const label* _ownStart,
            const label* _neiStart,
            const label* _own,
            const label* _nei,
            const label* _losort
        ):
             zero(_zero),
             Sf(_Sf),
             issf(_issf),
             ownStart(_ownStart),
             neiStart(_neiStart),
             own(_own),
             nei(_nei),
             losort(_losort)
        {}

        __HOST____DEVICE__
        GradType operator()(const label& id)
        {
            GradType out = zero;
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;

            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                out += Sf[face]*issf[face];
            }

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out -= Sf[face]*issf[face];
            }

            return out;
        }
    };

    template<class Type,class GradType>
    struct gaussGradPatchFunctor
    {
        const vector* Sf;
        const Type* issf;
        const label* neiStart;
        const label* losort;

        gaussGradPatchFunctor
        (
            const vector* _Sf,
            const Type* _issf,
            const label* _neiStart,
            const label* _losort
        ):
             Sf(_Sf),
             issf(_issf),
             neiStart(_neiStart),
             losort(_losort)
        {}

        __HOST____DEVICE__
        GradType operator()(const label& id,const GradType& g)
        {
            GradType out = g;
            
            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out += Sf[face]*issf[face];
            }

            return out;
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
Foam::fv::gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    const labelgpuList& l = mesh.lduAddr().lowerAddr();
    const labelgpuList& u = mesh.lduAddr().upperAddr();
    const labelgpuList& losort = mesh.lduAddr().losortAddr();

    const labelgpuList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = mesh.lduAddr().losortStartAddr();

    const vectorgpuField& Sf = mesh.Sf().getField();

    gpuField<GradType>& igGrad = gGrad.getField();
    const gpuField<Type>& issf = ssf.getField();
    
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+igGrad.size(),
        igGrad.begin(),
        gaussGradFunctor<Type,GradType>
        (
            pTraits<GradType>::zero,
            Sf.data(),
            issf.data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

    forAll(mesh.boundary(), patchi)
    {
        const vectorgpuField& pSf = mesh.Sf().boundaryField()[patchi];
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        
        const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchi);
        const labelgpuList& plosort = mesh.lduAddr().patchSortAddr(patchi);
        const labelgpuList& plosortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+pcells.size(),
            thrust::make_permutation_iterator
            (
                igGrad.begin(),
                pcells.begin()
            ),
            thrust::make_permutation_iterator(igGrad.begin(),pcells.begin()),
            gaussGradPatchFunctor<Type,GradType>
            (
                pSf.data(),
                pssf.data(),
                plosortStart.data(),
                plosort.data()
            )
        );

    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
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
Foam::fv::gaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        gradf(tinterpScheme_().interpolate(vsf), name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void Foam::fv::gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            const vectorgpuField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGrad.boundaryField()[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGrad.boundaryField()[patchi])
            );
        }
     }
}


// ************************************************************************* //
