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
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "pointConstraints.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

struct volPointInterpolationCalcPatchPointFunctor
{
    const label* faceNodes;
    bool* isPatchPoint;

    volPointInterpolationCalcPatchPointFunctor
    (
        const label* _faceNodes,
        bool* _isPatchPoint
    ):
        faceNodes(_faceNodes),
        isPatchPoint(_isPatchPoint)
    {}

    __host__ __device__
    void operator()(const faceData& face)
    {
        for(label i = face.start(); i < face.start()+face.size(); i++)
        {
            isPatchPoint[faceNodes[i]] = true;
        }
    } 
};

void volPointInterpolation::calcBoundaryAddressing()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::calcBoundaryAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    boundaryPtr_.reset
    (
        new primitivePatch
        (
            SubList<face>
            (
                mesh().faces(),
                mesh().nFaces()-mesh().nInternalFaces(),
                mesh().nInternalFaces()
            ),
            mesh().points()
        )
    );
    const primitivePatch& boundary = boundaryPtr_();

    boundaryIsPatchFace_.setSize(boundary.size());
    boundaryIsPatchFace_ = false;

    isPatchPoint_.setSize(mesh().nPoints());
    isPatchPoint_ = false;

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    // Get precalculated volField only so we can use coupled() tests for
    // cyclicAMI
    const surfaceScalarField& magSf = mesh().magSf();

    const faceDatagpuList& pFaces = boundary.getFaces();
    const labelgpuList& faceNodes = boundary.getFaceNodes();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if
        (
            !isA<emptyPolyPatch>(pp)
         && !magSf.boundaryField()[patchI].coupled()
        )
        {
            label pStart = pp.start()-mesh().nInternalFaces();
            label pEnd = pStart + pp.size();
            
            thrust::fill
            (
                boundaryIsPatchFace_.begin()+pStart,
                boundaryIsPatchFace_.begin()+pEnd,
                true
            );

            thrust::for_each
            (
                pFaces.begin()+pStart,
                pFaces.begin()+pEnd,
                volPointInterpolationCalcPatchPointFunctor
                (
                    faceNodes.data(),
                    isPatchPoint_.data()
                )
            );
        }
    }
}

struct volPointInterpolationMakeInternalWeightsFunctor
{
    const label* pointCellStart;
    const label* pointCells;
    const point* points;
    const vector* cellCentres;
    const bool* isPatchPoint;

    scalar* pointWeights;
    scalar* sumWeights;

    volPointInterpolationMakeInternalWeightsFunctor
    (
        const label* _pointCellStart,
        const label* _pointCells,
        const point* _points,
        const vector* _cellCentres,
        const bool* _isPatchPoint,

        scalar* _pointWeights,
        scalar* _sumWeights
    ):
        pointCellStart(_pointCellStart),
        pointCells(_pointCells),
        points(_points),
        cellCentres(_cellCentres),
        isPatchPoint(_isPatchPoint),
        pointWeights(_pointWeights),
        sumWeights(_sumWeights)
    {}

    __host__ __device__
    void operator()(const label& pointI)
    {
        if(!isPatchPoint[pointI])
        {
            scalar sum = 0;
            label start = pointCellStart[pointI];
            label end = pointCellStart[pointI+1];

            for(label i = start; i < end; i++)
            {
                scalar w = 1.0/mag(points[pointI] - cellCentres[pointCells[i]]);

                pointWeights[i] = w;
                sum += w;
            }

            sumWeights[pointI] = sum;
        }
    }
};

void volPointInterpolation::makeInternalWeights(scalargpuField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeInternalWeights() : "
            << "constructing weighting factors for internal and non-coupled"
            << " points." << endl;
    }

    const pointgpuField& points = mesh().getPoints();
    const labelgpuList& pointCells = mesh().getPointCells();
    const labelgpuList& pointCellsStart = mesh().getPointCellsStart();
    const vectorgpuField& cellCentres = mesh().getCellCentres();

    // Allocate storage for weighting factors
    pointWeights_.setSize(pointCells.size());

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    thrust::for_each
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(points.size()),
        volPointInterpolationMakeInternalWeightsFunctor
        (    
            pointCellsStart.data(),
            pointCells.data(),
            points.data(),
            cellCentres.data(),
            isPatchPoint_.data(),
            pointWeights_.data(),
            sumWeights.data()
        )
    );
}

struct volPointInterpolationMakeBoundaryWeightsFunctor
{
    const label nInternalFaces;
    const label* boundaryPoints;
    const label* pointFaceStart;
    const label* pointFaces;
    const point* points;
    const vector* faceCentres;
    const bool* isPatchPoint;
    const bool* boundaryIsPatchFace;

    scalar* pointWeights;
    scalar* sumWeights;

    volPointInterpolationMakeBoundaryWeightsFunctor
    (
        const label _nInternalFaces,
        const label* _boundaryPoints,
        const label* _pointFaceStart,
        const label* _pointFaces,
        const point* _points,
        const vector* _faceCentres,
        const bool* _isPatchPoint,
        const bool* _boundaryIsPatchFace,

        scalar* _pointWeights,
        scalar* _sumWeights
    ):
        nInternalFaces(_nInternalFaces),
        boundaryPoints(_boundaryPoints),
        pointFaceStart(_pointFaceStart),
        pointFaces(_pointFaces),
        points(_points),
        faceCentres(_faceCentres),
        isPatchPoint(_isPatchPoint),
        boundaryIsPatchFace(_boundaryIsPatchFace),
        pointWeights(_pointWeights),
        sumWeights(_sumWeights)
    {}

    __host__ __device__
    void operator()(const label& id)
    {
        label pointI = boundaryPoints[id];

        if (isPatchPoint[pointI])
        {
            scalar sum = 0;
            label start = pointFaceStart[id];
            label end = pointFaceStart[id+1];

            for(label i = start; i < end; i++)
            {
                label pFace = pointFaces[i];
                if (boundaryIsPatchFace[pFace])
                {
                    scalar w = 
                           1.0/mag(points[pointI] - faceCentres[nInternalFaces + pFace]);

                    pointWeights[i] = w;
                    sum += w;
                }
                else
                {
                    pointWeights[i] = 0.0;
                }
            }

            sumWeights[pointI] = sum;
        }
    }
};

void volPointInterpolation::makeBoundaryWeights(scalargpuField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeBoundaryWeights() : "
            << "constructing weighting factors for boundary points." << endl;
    }

    const pointgpuField& points = mesh().getPoints();
    const pointgpuField& faceCentres = mesh().getFaceCentres();

    const primitivePatch& boundary = boundaryPtr_();

    const labelgpuList& boundaryPoints = boundary.getMeshPoints();
    const labelgpuList& pointFaces = boundary.getPointFaces();
    const labelgpuList& pointFacesStart = boundary.getPointFacesStart();

    boundaryPointWeights_.setSize(pointFaces.size());

    thrust::for_each
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(boundaryPoints.size()),
        volPointInterpolationMakeBoundaryWeightsFunctor
        (
            mesh().nInternalFaces(),
            boundaryPoints.data(),
            pointFacesStart.data(),
            pointFaces.data(),
            points.data(),
            faceCentres.data(),
            isPatchPoint_.data(),
            boundaryIsPatchFace_.data(),
            boundaryPointWeights_.data(),
            sumWeights.data()
        )
    );
}

struct volPointInterpolationNormalizeWeights
{
    const label* pointDataStart;
    scalar * weights;

    volPointInterpolationNormalizeWeights
    (
        const label* _pointDataStart,
        scalar * _weights
    ):
        pointDataStart(_pointDataStart),
        weights(_weights)
    {}

    template<class Tuple>
    __host__ __device__
    void operator()(const Tuple& t)
    {
        const label& id = thrust::get<0>(t);
        const scalar sumWeight = thrust::get<1>(t);

        label start = pointDataStart[id];
        label end = pointDataStart[id+1];

        for(label i = start; i < end; i++)
        {
            weights[i] /= sumWeight;
        }
    }
};

void volPointInterpolation::makeWeights()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    // Update addressing over all boundary faces
    calcBoundaryAddressing();


    // Running sum of weights
    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar("zero", dimless, 0)
    );

    scalargpuField& sWeights = sumWeights.getField();


    // Create internal weights; add to sumWeights
    makeInternalWeights(sWeights);


    // Create boundary weights; override sumWeights
    makeBoundaryWeights(sWeights);


    //forAll(boundary.meshPoints(), i)
    //{
    //    label pointI = boundary.meshPoints()[i];
    //
    //    if (isPatchPoint_[pointI])
    //    {
    //        Pout<< "Calculated Weight at boundary point:" << i
    //            << " at:" << mesh().points()[pointI]
    //            << " sumWeight:" << sumWeights[pointI]
    //            << " from:" << boundaryPointWeights_[i]
    //            << endl;
    //    }
    //}


    // Sum collocated contributions
    pointConstraints::syncUntransformedData
    (
        mesh(),
        sWeights,
        plusEqOp<scalar>()
    );

    // And add separated contributions
    addSeparated(sumWeights);

    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves. Reuse the syncPointData
    // structure.
    pushUntransformedData(sWeights);


    // Normalise internal weights

    const labelgpuList& pointCellsStart = mesh().getPointCellsStart();
    thrust::for_each
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_counting_iterator(0),
            sWeights.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_counting_iterator(sumWeights.size()),
            sWeights.end()
        )),
        volPointInterpolationNormalizeWeights
        (
            pointCellsStart.data(),
            pointWeights_.data()
        )
    );

    // Normalise boundary weights
    const primitivePatch& boundary = boundaryPtr_();

    const labelgpuList& boundaryPoints = boundary.getMeshPoints();
    const labelgpuList& pointFacesStart = boundary.getPointFacesStart();

    thrust::for_each
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_counting_iterator(0),
            thrust::make_permutation_iterator
            (
                sWeights.begin(),
                boundaryPoints.begin()
            )
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_counting_iterator(boundaryPoints.size()),
            thrust::make_permutation_iterator
            (
                sWeights.begin(),
                boundaryPoints.end()
            )
        )),
        volPointInterpolationNormalizeWeights
        (
            pointFacesStart.data(),
            boundaryPointWeights_.data()
        )
    );

    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, volPointInterpolation>(vm)
{
    makeWeights();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

volPointInterpolation::~volPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volPointInterpolation::updateMesh(const mapPolyMesh&)
{
    makeWeights();
}


bool volPointInterpolation::movePoints()
{
    makeWeights();

    return true;
}


void volPointInterpolation::interpolateDisplacement
(
    const volVectorField& vf,
    pointVectorField& pf
) const
{
    interpolateInternalField(vf, pf);

    // Interpolate to the patches but no constraints
    interpolateBoundaryField(vf, pf);

    // Apply displacement constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrainDisplacement(pf, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
