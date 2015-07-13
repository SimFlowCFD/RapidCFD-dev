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

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "coupledPointPatchField.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void volPointInterpolation::pushUntransformedData
(
    gpuList<Type>& pointData
) const
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh().globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelgpuList& meshPoints = cpp.getMeshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    gpuList<Type> elemsGpu(elems.size());

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            pointData.begin(),
            meshPoints.begin()
        ),
        thrust::make_permutation_iterator
        (
            pointData.begin(),
            meshPoints.end()
        ),
        elemsGpu.begin()
    );

    thrust::copy
    (
        elemsGpu.begin(),
        elemsGpu.end(),
        elems.begin()
    );

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    elemsGpu = elems;
    thrust::copy
    (
        elemsGpu.begin(),
        elemsGpu.end(),
        thrust::make_permutation_iterator
        (
            pointData.begin(),
            meshPoints.begin()
        )
    );
}


template<class Type>
void volPointInterpolation::addSeparated
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::addSeparated" << endl;
    }

    forAll(pf.boundaryField(), patchI)
    {
        if (pf.boundaryField()[patchI].coupled())
        {
            refCast<coupledPointPatchField<Type> >
                (pf.boundaryField()[patchI]).initSwapAddSeparated
                (
                    Pstream::nonBlocking,
                    pf.internalField()
                );
        }
    }

    // Block for any outstanding requests
    Pstream::waitRequests();

    forAll(pf.boundaryField(), patchI)
    {
        if (pf.boundaryField()[patchI].coupled())
        {
            refCast<coupledPointPatchField<Type> >
                (pf.boundaryField()[patchI]).swapAddSeparated
                (
                    Pstream::nonBlocking,
                    pf.internalField()
                );
        }
    }
}

template<class Type>
struct volInterpolationInterpolationInternalFieldFunctor
{
    const label* pointCellsStart;
    const label* pointCells;
    const scalar* weight;
    const bool* isPatchPoint;
    const Type* vf;
    Type* pf;
    const Type zero;

    volInterpolationInterpolationInternalFieldFunctor
    (
        const label* _pointCellsStart,
        const label* _pointCells,
        const scalar* _weight,
        const bool* _isPatchPoint,
        const Type* _vf,
        Type* _pf
    ):
        pointCellsStart(_pointCellsStart),
        pointCells(_pointCells),
        weight(_weight),
        isPatchPoint(_isPatchPoint),
        vf(_vf),
        pf(_pf),
        zero(pTraits<Type>::zero)
    {}

    __host__ __device__
    void operator()(const label& pointI)
    {
        if(!isPatchPoint[pointI])
        {
            const label start = pointCellsStart[pointI];
            const label end = pointCellsStart[pointI+1];

            Type out = zero;

            for(label i = start; i < end; i++)
            {
                out += weight[i]*vf[pointCells[i]];
            }
 
            pf[pointI] = out;
        }
    }
};

template<class Type>
void volPointInterpolation::interpolateInternalField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateInternalField("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelgpuList& pointCells = vf.mesh().getPointCells();
    const labelgpuList& pointCellsStart = vf.mesh().getPointCellsStart();

    // Multiply volField by weighting factor matrix to create pointField
    thrust::for_each
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(pf.size()),
        volInterpolationInterpolationInternalFieldFunctor<Type>
        (
            pointCellsStart.data(),
            pointCells.data(),
            pointWeights_.data(),
            isPatchPoint_.data(),
            vf.getField().data(),
            pf.getField().data()
        )
    );
}


template<class Type>
tmp<gpuField<Type> > volPointInterpolation::flatBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const fvBoundaryMesh& bm = mesh.boundary();

    tmp<gpuField<Type> > tboundaryVals
    (
        new gpuField<Type>(mesh.nFaces()-mesh.nInternalFaces())
    );
    gpuField<Type>& boundaryVals = tboundaryVals();

    forAll(vf.boundaryField(), patchI)
    {
        const polyPatch& pp = bm[patchI].patch();
        label bFaceI = pp.start() - mesh.nInternalFaces();

        if
        (
           !isA<emptyFvPatch>(pp)
        && !vf.boundaryField()[patchI].coupled()
        )
        {
            gpuList<Type>
            (
                boundaryVals,
                vf.boundaryField()[patchI].size(),
                bFaceI
            ) = vf.boundaryField()[patchI];
        }
        else
        {
            thrust::fill
            (
                boundaryVals.begin()+bFaceI,
                boundaryVals.begin()+bFaceI+pp.size(),
                pTraits<Type>::zero
            );
        }
    }

    return tboundaryVals;
}

template<class Type>
struct volPointInterpolationInterpolateBoundaryFieldFunctor
{
    const label* points;
    const label* pointFacesStart;
    const label* pointFaces;
    const bool* isPatchPoint;
    const bool* boundaryIsPatchFace;
    const scalar* weight;
    const Type* bf;
    Type* pf;
    const Type zero;
    
    volPointInterpolationInterpolateBoundaryFieldFunctor
    (
        const label* _points,
        const label* _pointFacesStart,
        const label* _pointFaces,
        const bool* _isPatchPoint,
        const bool* _boundaryIsPatchFace,
        const scalar* _weight,
        const Type* _bf,
        Type* _pf
    ):
        points(_points),
        pointFacesStart(_pointFacesStart),
        pointFaces(_pointFaces),
        isPatchPoint(_isPatchPoint),
        boundaryIsPatchFace(_boundaryIsPatchFace),
        weight(_weight),
        bf(_bf),
        pf(_pf),
        zero(pTraits<Type>::zero)
    {}

    __host__ __device__
    void operator()(const label& id)
    {
        const label pointI = points[id];

        if (isPatchPoint[pointI])
        {
            const label start = pointFacesStart[id];
            const label end = pointFacesStart[id+1];

            Type out = zero;

            for(label i = start; i < end; i++)
            {
                label faceI = pointFaces[i];

                if (boundaryIsPatchFace[faceI])
                {
                   out += weight[i]*bf[faceI];
                }
            }

            pf[pointI] = out;
        }
    }
};

template<class Type>
void volPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    const primitivePatch& boundary = boundaryPtr_();

    gpuField<Type>& pfi = pf.internalField();

    // Get face data in flat list
    tmp<gpuField<Type> > tboundaryVals(flatBoundaryField(vf));
    const gpuField<Type>& boundaryVals = tboundaryVals();

    const labelgpuList& boundaryPoints = boundary.getMeshPoints();
    const labelgpuList& pointFaces = boundary.getPointFaces();
    const labelgpuList& pointFacesStart = boundary.getPointFacesStart();

    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    thrust::for_each
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(boundaryPoints.size()),
        volPointInterpolationInterpolateBoundaryFieldFunctor<Type>
        (
            boundaryPoints.data(),
            pointFacesStart.data(),
            pointFaces.data(),
            isPatchPoint_.data(),
            boundaryIsPatchFace_.data(),
            boundaryPointWeights_.data(),
            boundaryVals.data(),
            pfi.data()
        )
    );

    // Sum collocated contributions
    pointConstraints::syncUntransformedData(mesh(), pfi, plusEqOp<Type>());

    // And add separated contributions
    addSeparated(pf);

    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves.
    pushUntransformedData(pfi);
}


template<class Type>
void volPointInterpolation::interpolateBoundaryField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const bool overrideFixedValue
) const
{
    interpolateBoundaryField(vf, pf);

    // Apply constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrain(pf, overrideFixedValue);
}


template<class Type>
void volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    interpolateInternalField(vf, pf);

    // Interpolate to the patches preserving fixed value BCs
    interpolateBoundaryField(vf, pf, false);
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const wordList& patchFieldTypes
) const
{
    const pointMesh& pm = pointMesh::New(vf.mesh());

    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "volPointInterpolate(" + vf.name() + ')',
                vf.instance(),
                pm.thisDb()
            ),
            pm,
            vf.dimensions(),
            patchFieldTypes
        )
    );

    interpolateInternalField(vf, tpf());

    // Interpolate to the patches overriding fixed value BCs
    interpolateBoundaryField(vf, tpf(), true);

    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const wordList& patchFieldTypes
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf(), patchFieldTypes);
    tvf.clear();
    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name,
    const bool cache
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    const pointMesh& pm = pointMesh::New(vf.mesh());
    const objectRegistry& db = pm.thisDb();

    if (!cache || vf.mesh().changing())
    {
        // Delete any old occurences to avoid double registration
        if (db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            PointFieldType& pf = const_cast<PointFieldType&>
            (
                db.objectRegistry::template lookupObject<PointFieldType>(name)
            );

            if (pf.ownedByRegistry())
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;
            }
        }


        tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
        (
            new GeometricField<Type, pointPatchField, pointMesh>
            (
                IOobject
                (
                    name,
                    vf.instance(),
                    pm.thisDb()
                ),
                pm,
                vf.dimensions()
            )
        );

        interpolate(vf, tpf());

        return tpf;
    }
    else
    {
        if (!db.objectRegistry::template foundObject<PointFieldType>(name))
        {
            solution::cachePrintMessage("Calculating and caching", name, vf);
            tmp<PointFieldType> tpf = interpolate(vf, name, false);
            PointFieldType* pfPtr = tpf.ptr();
            regIOobject::store(pfPtr);
            return *pfPtr;
        }
        else
        {
            PointFieldType& pf = const_cast<PointFieldType&>
            (
                db.objectRegistry::template lookupObject<PointFieldType>(name)
            );

            if (pf.upToDate(vf))    //TBD: , vf.mesh().points()))
            {
                solution::cachePrintMessage("Reusing", name, vf);
                return pf;
            }
            else
            {
                solution::cachePrintMessage("Deleting", name, vf);
                pf.release();
                delete &pf;

                solution::cachePrintMessage("Recalculating", name, vf);
                tmp<PointFieldType> tpf = interpolate(vf, name, false);

                solution::cachePrintMessage("Storing", name, vf);
                PointFieldType* pfPtr = tpf.ptr();
                regIOobject::store(pfPtr);

                // Note: return reference, not pointer
                return *pfPtr;
            }
        }
    }
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return interpolate(vf, "volPointInterpolate(" + vf.name() + ')', false);
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf());
    tvf.clear();
    return tpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
