/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "directAMI.H"
#include "primitiveFieldsFwd.H"
#include "pointFieldFwd.H"
#include "faceData.H"
#include "faceFunctors.H"
#include "gpuList.H"
#include "scalarList.H"

namespace Foam {
    struct ShootRayFunctor {
        const point *tgtPoints;
        const point *srcPoints;
        const vector *srcCf;
        const faceData *tgtFaces;
        const faceData *srcFaces;
        const label size;

        // Arrays used for addressing
        label* srcAddr;
        label* tgtAddr;

        // Arrays used to hold the actual count
        label* srcCount;
        label* tgtCount;

        faceRayFunctor rayFunctor;
        faceNormalFunctor normalFunctor;

        ShootRayFunctor
        (
            const point *_tgtPoints,
            const point *_srcPoints,
            const vector *_srcCf,
            const faceData *_tgtFaces,
            const faceData *_srcFaces,
            const label _size,
            label* _srcCount,
            label* _tgtCount,
            label* _srcAddr,
            label* _tgtAddr
        ) :
            tgtPoints(_tgtPoints),
            srcPoints(_srcPoints),
            srcCf(_srcCf),
            tgtFaces(_tgtFaces),
            srcFaces(_srcFaces),
            size(_size),
            srcCount(_srcCount),
            tgtCount(_tgtCount),
            srcAddr(_srcAddr),
            tgtAddr(_tgtAddr),
            rayFunctor(_tgtPoints, intersection::FULL_RAY, intersection::VECTOR),
            normalFunctor(_srcPoints)
        {}

        template <typename Tuple>
        __host__ __device__
        void operator()(Tuple t)
        {
            const label tgtI = thrust::get<0>(t);
            const label srcI = thrust::get<1>(t);

            if (rayFunctor(tgtFaces[tgtI], srcCf[srcI], normalFunctor(srcFaces[srcI])).hit()) {
                if (srcCount[srcI] < size)
                {
                    srcAddr[srcI*size + srcCount[srcI]] = tgtI;
                    srcCount[srcI]++;
                }
                if (tgtCount[tgtI] < size)
                {
                    tgtAddr[tgtI*size + tgtCount[tgtI]] = srcI;
                    tgtCount[tgtI]++;
                }
            }
        }
    };
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::appendToDirectSeeds
(
    labelList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    DynamicList<label>& nonOverlapFaces,
    label& srcFaceI,
    label& tgtFaceI
) const
{
    const labelList& srcNbr = this->srcPatch_.faceFaces()[srcFaceI];
    const labelList& tgtNbr = this->tgtPatch_.faceFaces()[tgtFaceI];

    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    const vectorField& srcCf = this->srcPatch_.faceCentres();

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if ((mapFlag[srcI] == 0) && (srcTgtSeed[srcI] == -1))
        {
            // first attempt: match by comparing face centres
            const face& srcF = this->srcPatch_[srcI];
            const point& srcC = srcCf[srcI];

            scalar tol = GREAT;
            forAll(srcF, fpI)
            {
                const point& p = srcPoints[srcF[fpI]];
                scalar d2 = magSqr(p - srcC);
                if (d2 < tol)
                {
                    tol = d2;
                }
            }
            tol = max(SMALL, 0.0001*sqrt(tol));

            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];
                const face& tgtF = this->tgtPatch_[tgtI];
                const point tgtC = tgtF.centre(tgtPoints);

                if (mag(srcC - tgtC) < tol)
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            // second attempt: match by shooting a ray into the tgt face
            if (!found)
            {
                const vector srcN = srcF.normal(srcPoints);

                forAll(tgtNbr, j)
                {
                    label tgtI = tgtNbr[j];
                    const face& tgtF = this->tgtPatch_[tgtI];
                    pointHit ray = tgtF.ray(srcCf[srcI], srcN, tgtPoints);

                    if (ray.hit())
                    {
                        // new match - append to lists
                        found = true;

                        srcTgtSeed[srcI] = tgtI;
                        srcSeeds.append(srcI);

                        break;
                    }
                }
            }

            // no match available for source face srcI
            if (!found)
            {
                mapFlag[srcI] = -1;
                nonOverlapFaces.append(srcI);

                if (debug)
                {
                    Pout<< "source face not found: id=" << srcI
                        << " centre=" << srcCf[srcI]
                        << " normal=" << srcF.normal(srcPoints)
                        << " points=" << srcF.points(srcPoints)
                        << endl;

                    Pout<< "target neighbours:" << nl;
                    forAll(tgtNbr, j)
                    {
                        label tgtI = tgtNbr[j];
                        const face& tgtF = this->tgtPatch_[tgtI];

                        Pout<< "face id: " << tgtI
                            << " centre=" << tgtF.centre(tgtPoints)
                            << " normal=" << tgtF.normal(tgtPoints)
                            << " points=" << tgtF.points(tgtPoints)
                            << endl;
                    }
                }
            }
        }
    }

    if (srcSeeds.size())
    {
        srcFaceI = srcSeeds.remove();
        tgtFaceI = srcTgtSeed[srcFaceI];
    }
    else
    {
        srcFaceI = -1;
        tgtFaceI = -1;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::restartAdvancingFront
(
    labelList& mapFlag,
    DynamicList<label>& nonOverlapFaces,
    label& srcFaceI,
    label& tgtFaceI
) const
{
    forAll(mapFlag, faceI)
    {
        if (mapFlag[faceI] == 0)
        {
            tgtFaceI = this->findTargetFace(faceI);

            if (tgtFaceI < 0)
            {
                mapFlag[faceI] = -1;
                nonOverlapFaces.append(faceI);
            }
            else
            {
                srcFaceI = faceI;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::directAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::~directAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFaceI,
    label tgtFaceI
)
{
    this->checkPatches();

    // set initial sizes for weights and addressing
    srcAddress.setSize(this->srcPatch_.size());
    srcWeights.setSize(this->srcPatch_.size());
    tgtAddress.setSize(this->tgtPatch_.size());
    tgtWeights.setSize(this->tgtPatch_.size());

    // check that patch sizes are valid
    if (!this->srcPatch_.size())
    {
        return;
    }
    else if (!this->tgtPatch_.size())
    {
        WarningIn
            (
                "void Foam::AMIMethod<SourcePatch, TargetPatch>::initialise"
                    "("
                    "labelListList&, "
                    "scalarListList&, "
                    "labelListList&, "
                    "scalarListList&, "
                    "label&, "
                    "label&"
                    ")"
            )
        << this->srcPatch_.size() << " source faces but no target faces" << endl;

        return;
    }

    const pointgpuField& tgtPoints = this->tgtPatch_.getPoints();
    const pointgpuField& srcPoints = this->srcPatch_.getPoints();

    const faceDatagpuList& tgtFaces = this->tgtPatch_.getFaces();
    const faceDatagpuList& srcFaces = this->srcPatch_.getFaces();

    vectorgpuField srcCf(srcFaces.size());

    thrust::transform
    (
        srcFaces.begin(),
        srcFaces.end(),
        srcCf.begin(),
        faceCentreFunctor
        (
            this->srcPatch_.getFaceNodes().data(),
            srcPoints.data()
        )
    );

    label size = 15;

    labelgpuList srcCount(this->srcPatch_.size(), 0);
    labelgpuList tgtCount(this->tgtPatch_.size(), 0);

    labelgpuList srcAddr(this->srcPatch_.size()*size, -1);
    labelgpuList tgtAddr(this->tgtPatch_.size()*size, -1);

    ShootRayFunctor shootRayFunctor
        (
            tgtPoints.data(),
            srcPoints.data(),
            srcCf.data(),
            tgtFaces.data(),
            srcFaces.data(),
            size,
            srcCount.data(),
            tgtCount.data(),
            srcAddr.data(),
            tgtAddr.data()
        );

    for (label i = 0; i < srcFaces.size(); ++i)
    {
        thrust::for_each
        (
            thrust::make_zip_iterator
            (
                thrust::make_tuple
                (
                    thrust::counting_iterator<label>(0),
                    thrust::constant_iterator<label>(i)
                )
            ),
            thrust::make_zip_iterator
            (
                thrust::make_tuple
                (
                    thrust::counting_iterator<label>(tgtFaces.size()),
                    thrust::constant_iterator<label>(i)
                )
            ),
            shootRayFunctor
        );
    }

    for (label i = 0; i < this->srcPatch_.size(); ++i)
    {
        scalar magSf = this->srcMagSf_[i];
        srcWeights[i] = scalarList(1, magSf);

        labelList l(srcCount.get(i));

        for (label j = 0; j < srcCount.get(i); ++j)
        {
            l[j] = srcAddr.get(size*i + j);
        }

        srcAddress[i].transfer(l);
    }

    for (label i = 0; i < this->tgtPatch_.size(); ++i)
    {
        scalar magSf = this->tgtMagSf_[i];
        tgtWeights[i] = scalarList(1, magSf);

        labelList l(tgtCount.get(i));

        for (label j = 0; j < tgtCount.get(i); ++j)
        {
            l[j] = tgtAddr.get(size*i + j);
        }

        tgtAddress[i].transfer(l);
    }
}
// ************************************************************************* //
