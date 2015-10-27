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

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

    template<class Type,bool integrate>
    struct surfaceIntegrateFunctor : public std::unary_function<label,Type>
    {
        const Type zero;
        const Type* issf;
        const label* ownStart;
        const label* neiStart;
        const label* own;
        const label* nei;
        const label* losort;

        surfaceIntegrateFunctor
        (
             const Type* _issf,
             const label* _ownStart,
             const label* _neiStart,
             const label* _own,
             const label* _nei,
             const label* _losort
        ):
             zero(pTraits<Type>::zero),
             issf(_issf),
             ownStart(_ownStart),
             neiStart(_neiStart),
             own(_own),
             nei(_nei),
             losort(_losort)
        {}

        __HOST____DEVICE__
        Type operator()(const label& id)
        {
            Type out = zero;
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;

            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                out += issf[face];
            }

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                if(integrate)
                    out -= issf[face];
                else
                   out += issf[face];
            }

            return out;
        }
    };

    template<class Type>
    struct surfaceIntegratePatchFunctor : public std::binary_function<label,Type,Type>
    {
        const Type* issf;
        const label* neiStart;
        const label* losort;

        surfaceIntegratePatchFunctor
        (
            const Type* _issf,
            const label* _neiStart,
            const label* _losort
        ):
             issf(_issf),
             neiStart(_neiStart),
             losort(_losort)
        {}

        __HOST____DEVICE__
        Type operator()(const label& id,const Type& t)
        {
            Type out = t;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out += issf[face];
            }

            return out;
        }
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    gpuField<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelgpuList& l = mesh.lduAddr().lowerAddr();
    const labelgpuList& u = mesh.lduAddr().upperAddr();
    const labelgpuList& losort = mesh.lduAddr().losortAddr();

    const labelgpuList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = mesh.lduAddr().losortStartAddr();

    const gpuField<Type>& issf = ssf.getField();

    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+ivf.size(),
        ivf.begin(),
        surfaceIntegrateFunctor<Type,true>
        (
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
                ivf.begin(),
                pcells.begin()
            ),
            thrust::make_permutation_iterator
            (
                ivf.begin(),
                pcells.begin()
            ),
            surfaceIntegratePatchFunctor<Type>
            (
                pssf.data(),
                plosortStart.data(),
                plosort.data()
            )
        );
    }

    ivf /= mesh.Vsc()().getField();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimVol,
                pTraits<Type>::zero
            ),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    surfaceIntegrate(vf.internalField(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    const labelgpuList& l = mesh.lduAddr().lowerAddr();
    const labelgpuList& u = mesh.lduAddr().upperAddr();
    const labelgpuList& losort = mesh.lduAddr().losortAddr();

    const labelgpuList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = mesh.lduAddr().losortStartAddr();


    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+vf.size(),
        vf.getField().begin(),
        surfaceIntegrateFunctor<Type,false>
        (
            ssf.getField().data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

/*
    forAll(owner, facei)
    {
        vf[owner[facei]] += ssf[facei];
        vf[neighbour[facei]] += ssf[facei];
    }
*/
    forAll(mesh.boundary(), patchi)
    {
        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchi);

        const labelgpuList& losort = mesh.lduAddr().patchSortAddr(patchi);
        const labelgpuList& losortStart = mesh.lduAddr().patchSortStartAddr(patchi);
/*
        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += pssf[facei];
        }
*/

        thrust::transform
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+pcells.size(),
            thrust::make_permutation_iterator
            (
                vf.getField().begin(),
                pcells.begin()
            ),
            thrust::make_permutation_iterator
            (
                vf.getField().begin(),
                pcells.begin()
            ),
            surfaceIntegratePatchFunctor<Type>
            (
                pssf.data(),
                losortStart.data(),
                losort.data()
            )
        );
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
