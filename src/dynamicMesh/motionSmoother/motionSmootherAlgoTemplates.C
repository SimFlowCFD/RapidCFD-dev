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

#include "motionSmootherAlgo.H"
#include "meshTools.H"
#include "processorPointPatchFields.H"
#include "pointConstraint.H"
#include "pointConstraints.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
void Foam::motionSmootherAlgo::checkConstraints
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> FldType;

    const polyMesh& mesh = pf.mesh();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyPolyPatch>(bm[patchi]))
        {
            nPatchPatchPoints += bm[patchi].boundaryPoints().size();
        }
    }


    typename FldType::GeometricBoundaryField& bFld = pf.boundaryField();


    // Evaluate in reverse order

    forAllReverse(bFld, patchi)
    {
        bFld[patchi].initEvaluate(Pstream::blocking);   // buffered
    }

    forAllReverse(bFld, patchi)
    {
        bFld[patchi].evaluate(Pstream::blocking);
    }


    // Save the values

    Field<Type> boundaryPointValues(nPatchPatchPoints);
    nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyPolyPatch>(bm[patchi]))
        {
            const labelList& bp = bm[patchi].boundaryPoints();
            const labelList& meshPoints = bm[patchi].meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];
                boundaryPointValues[nPatchPatchPoints++] = pf[ppp];
            }
        }
    }


    // Forward evaluation

    bFld.evaluate();


    // Check

    nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if (!isA<emptyPolyPatch>(bm[patchi]))
        {
            const labelList& bp = bm[patchi].boundaryPoints();
            const labelList& meshPoints = bm[patchi].meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                const Type& savedVal = boundaryPointValues[nPatchPatchPoints++];

                if (savedVal != pf[ppp])
                {
                    FatalErrorIn
                    (
                        "motionSmootherAlgo::checkConstraints"
                        "(GeometricField<Type, pointPatchField, pointMesh>&)"
                    )   << "Patch fields are not consistent on mesh point "
                        << ppp << " coordinate " << mesh.points()[ppp]
                        << " at patch " << bm[patchi].name() << '.'
                        << endl
                        << "Reverse evaluation gives value " << savedVal
                        << " , forward evaluation gives value " << pf[ppp]
                        << abort(FatalError);
                }
            }
        }
    }
}

namespace Foam
{

struct motionSmootherAlgoEdgeNode
{
    const int index;

    motionSmootherAlgoEdgeNode
    (
        const int _index
    ):
        index(_index)
    {}

    __host__ __device__
    label operator()(const edge& e)
    {
        return e[index];
    }
};

template<class Type>
struct motionSmootherAlgoWeightsAvgFunctor
{
    const bool first;
    const __restrict__ Type * fld;
    const Type zero;

    motionSmootherAlgoWeightsAvgFunctor
    (
        const bool _first,
        const __restrict__ Type * _fld
    ):
        first(_first),
        fld(_fld),
        zero(pTraits<Type>::zero)
    {}
    
    template<class Tuple>
    __host__ __device__
    thrust::tuple<Type,scalar>
    operator()(const bool& isMaster, const Tuple& t)
    {
        if(isMaster)
        {
            const edge& e = thrust::get<0>(t);
            const scalar w = thrust::get<1>(t);

            return thrust::make_tuple
                   (
                       w*fld[e[first?1:0]],
                       w
                   );
        }
        else
        {
            return thrust::make_tuple(zero,0.0);
        }
    }
};

template<class Type>
struct motionSmootherAlgoAvgFunctor
{
    template<class Tuple>
    __host__ __device__
    Type operator()(const Type& res, const Tuple& t)
    {
        if (mag(thrust::get<0>(t)) < VSMALL)
        {
            // Unconnected point. Take over original value
            return thrust::get<1>(t);
        }
        else
        {
            return res/thrust::get<0>(t);
        }
    }
};

};

// Average of connected points.
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::motionSmootherAlgo::avg
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalargpuField& edgeWeight
) const
{
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tres
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "avg("+fld.name()+')',
                fld.time().timeName(),
                fld.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fld.mesh(),
            dimensioned<Type>("zero", fld.dimensions(), pTraits<Type>::zero)
        )
    );
    GeometricField<Type, pointPatchField, pointMesh>& res = tres();

    const polyMesh& mesh = fld.mesh()();


    // Sum local weighted values and weights
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note: on coupled edges use only one edge (through isMasterEdge)
    // This is done so coupled edges do not get counted double.

    scalargpuField sumWeight(mesh.nPoints(), 0.0);

    const edgegpuList& edges = mesh.getEdges();

    labelgpuList nodeTmp(edges.size());
    gpuList<Type> resTmp(edges.size());
    gpuList<scalar> sumWeightTmp(edges.size());

    labelgpuList nodeOut(edges.size());
    gpuList<Type> resOut(edges.size());
    scalargpuList sumWeightOut(edges.size());

    // begin edge
    thrust::transform
    (
        gpuIsMasterEdge_.begin(),
        gpuIsMasterEdge_.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            edges.begin(),
            edgeWeight.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            resTmp.begin(),
            sumWeightTmp.begin()
        )),
        motionSmootherAlgoWeightsAvgFunctor<Type>
        (
            true,
            fld.getField().data()
        )
    );

    // reduce res - begin
    thrust::transform
    (
        edges.begin(),
        edges.end(),
        nodeTmp.begin(),
        motionSmootherAlgoEdgeNode(0)
    );

    thrust::sort_by_key
    (
        nodeTmp.begin(),
        nodeTmp.end(),
        resTmp.begin()
    );

    typename thrust::pair
             <
                 typename gpuList<label>::iterator,
                 typename gpuList<Type>::iterator
             > resEnd = thrust::reduce_by_key
             (
                  nodeTmp.begin(),
                  nodeTmp.end(),
                  resTmp.begin(),
                  nodeOut.begin(),
                  resOut.begin()
             );

    thrust::copy
    (
        resOut.begin(),
        resEnd.second,
        thrust::make_permutation_iterator
        (
            res.getField().begin(),
            nodeOut.begin()
        )
    );

    // reduce sumWeights - begin
    thrust::transform
    (
        edges.begin(),
        edges.end(),
        nodeTmp.begin(),
        motionSmootherAlgoEdgeNode(0)
    );

    thrust::sort_by_key
    (
        nodeTmp.begin(),
        nodeTmp.end(),
        sumWeightTmp.begin()
    );

    typename thrust::pair
             <
                 typename gpuList<label>::iterator,
                 typename gpuList<scalar>::iterator
             > sumWeightEnd = thrust::reduce_by_key
             (
                  nodeTmp.begin(),
                  nodeTmp.end(),
                  sumWeightTmp.begin(),
                  nodeOut.begin(),
                  sumWeightOut.begin()
             );

    thrust::copy
    (
        sumWeightOut.begin(),
        sumWeightEnd.second,
        thrust::make_permutation_iterator
        (
            sumWeight.begin(),
            nodeOut.begin()
        )
    );

    // end edge

    sumWeightTmp = 0.0;
    sumWeightOut = 0.0;
    resTmp = pTraits<Type>::zero;
    resOut = pTraits<Type>::zero;
     
    thrust::transform
    (
        gpuIsMasterEdge_.begin(),
        gpuIsMasterEdge_.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            edges.begin(),
            edgeWeight.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            resTmp.begin(),
            sumWeightTmp.begin()
        )),
        motionSmootherAlgoWeightsAvgFunctor<Type>
        (
            false,
            fld.getField().data()
        )
    );

    // reduce res - end
    thrust::transform
    (
        edges.begin(),
        edges.end(),
        nodeTmp.begin(),
        motionSmootherAlgoEdgeNode(1)
    );

    thrust::sort_by_key
    (
        nodeTmp.begin(),
        nodeTmp.end(),
        resTmp.begin()
    );

    resEnd = thrust::reduce_by_key
             (
                  nodeTmp.begin(),
                  nodeTmp.end(),
                  resTmp.begin(),
                  nodeOut.begin(),
                  resOut.begin()
             );

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            res.getField().begin(),
            nodeOut.begin()
        ),
        thrust::make_permutation_iterator
        (
            res.getField().begin(),
            resEnd.first
        ),
        resOut.begin(),
        thrust::make_permutation_iterator
        (
            res.getField().begin(),
            nodeOut.begin()
        ),
        thrust::plus<Type>()
    );

    // reduce sumWeights - begin
    thrust::transform
    (
        edges.begin(),
        edges.end(),
        nodeTmp.begin(),
        motionSmootherAlgoEdgeNode(1)
    );

    thrust::sort_by_key
    (
        nodeTmp.begin(),
        nodeTmp.end(),
        sumWeightTmp.begin()
    );

    sumWeightEnd = thrust::reduce_by_key
             (
                  nodeTmp.begin(),
                  nodeTmp.end(),
                  sumWeightTmp.begin(),
                  nodeOut.begin(),
                  sumWeightOut.begin()
             );

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            sumWeight.begin(),
            nodeOut.begin()
        ),
        thrust::make_permutation_iterator
        (
            sumWeight.begin(),
            sumWeightEnd.first
        ),
        sumWeightOut.begin(),
        thrust::make_permutation_iterator
        (
            sumWeight.begin(),
            nodeOut.begin()
        ),
        thrust::plus<scalar>()
    );

/*
    forAll(edges, edgeI)
    {
        if (isMasterEdge_.get(edgeI) == 1)
        {
            const edge& e = edges[edgeI];
            const scalar w = edgeWeight[edgeI];

            res[e[0]] += w*fld[e[1]];
            sumWeight[e[0]] += w;

            res[e[1]] += w*fld[e[0]];
            sumWeight[e[1]] += w;
        }
    }
*/

    // Add coupled contributions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    syncTools::syncPointList
    (
        mesh,
        res,
        plusEqOp<Type>(),
        pTraits<Type>::zero     // null value
    );
    syncTools::syncPointList
    (
        mesh,
        sumWeight,
        plusEqOp<scalar>(),
        scalar(0)               // null value
    );


    // Average
    // ~~~~~~~

    thrust::transform
    (
        res.getField().begin(),
        res.getField().end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            sumWeight.begin(),
            fld.getField().begin()
        )),
        res.getField().begin(),
        motionSmootherAlgoAvgFunctor<Type>()
    );
/*    forAll(res, pointI)
    {
        if (mag(sumWeight[pointI]) < VSMALL)
        {
            // Unconnected point. Take over original value
            res[pointI] = fld[pointI];
        }
        else
        {
            res[pointI] /= sumWeight[pointI];
        }
    }*/

    // Single and multi-patch constraints
    pointConstraints::New(fld.mesh()).constrain(res, false);

    return tres;
}

namespace Foam
{

template<class Type>
struct motionSmootherAlgoSmoothFunctor
{
    template<class Tuple>
    point operator()(const bool& isInternal, const Tuple& t)
    {
        const point& newFld = thrust::get<0>(t);
        const point& fld = thrust::get<1>(t);
        const point& avgFld = thrust::get<2>(t);

        if(isInternal)
        {
            return 0.5*fld + 0.5*avgFld;
        }
        else
        {
            return newFld;
        }       
    }
};

}

// smooth field (point-jacobi)
template<class Type>
void Foam::motionSmootherAlgo::smooth
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalarField& edgeWeight,
    GeometricField<Type, pointPatchField, pointMesh>& newFld
) const
{
    tmp<pointVectorField> tavgFld = avg(fld, edgeWeight);
    const pointVectorField& avgFld = tavgFld();

    thrust::transform
    (
        gpuIsInternalPoint_.begin(),
        gpuIsInternalPoint_.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            newFld.getField().begin(),
            fld.getField().begin(),
            avgFld.getField().begin()
        )),
        newFld.getField().begin(), 
        motionSmootherAlgoSmoothFunctor<Type>()
    );
/*
    forAll(fld, pointI)
    {
        if (isInternalPoint(pointI))
        {
            newFld[pointI] = 0.5*fld[pointI] + 0.5*avgFld[pointI];
        }
    }
*/
    // Single and multi-patch constraints
    pointConstraints::New(fld.mesh()).constrain(newFld, false);
}


//- Test synchronisation of generic field (not positions!) on points
template<class Type, class CombineOp>
void Foam::motionSmootherAlgo::testSyncField
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const Type& zero,
    const scalar maxMag
) const
{
    if (debug)
    {
        Pout<< "testSyncField : testing synchronisation of Field<Type>."
            << endl;
    }

    Field<Type> syncedFld(fld);

    syncTools::syncPointList
    (
        mesh_,
        syncedFld,
        cop,
        zero
    );

    forAll(syncedFld, i)
    {
        if (mag(syncedFld[i] - fld[i]) > maxMag)
        {
            FatalErrorIn
            (
                "motionSmootherAlgo::testSyncField"
                "(const Field<Type>&, const CombineOp&"
                ", const Type&, const bool)"
            )   << "On element " << i << " value:" << fld[i]
                << " synchronised value:" << syncedFld[i]
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
