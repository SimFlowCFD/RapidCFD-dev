/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "Time.H"
#include "GAMGInterface.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
    defineRunTimeSelectionTable(GAMGAgglomeration, geometry);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels(const label nCreatedLevels)
{
    nCells_.setSize(nCreatedLevels);
    restrictSortAddressing_.setSize(nCreatedLevels),
    restrictTargetAddressing_.setSize(nCreatedLevels),
    restrictTargetStartAddressing_.setSize(nCreatedLevels),
    restrictAddressingHost_.setSize(nCreatedLevels),
    nFaces_.setSize(nCreatedLevels);
    faceRestrictSortAddressing_.setSize(nCreatedLevels),
    faceRestrictTargetAddressing_.setSize(nCreatedLevels),
    faceRestrictTargetStartAddressing_.setSize(nCreatedLevels),
    faceRestrictAddressingHost_.setSize(nCreatedLevels),
    faceFlipMap_.setSize(nCreatedLevels);
    faceFlipMapHost_.setSize(nCreatedLevels);
    nPatchFaces_.setSize(nCreatedLevels);
    patchFaceRestrictSortAddressing_.setSize(nCreatedLevels),
    patchFaceRestrictTargetAddressing_.setSize(nCreatedLevels),
    patchFaceRestrictTargetStartAddressing_.setSize(nCreatedLevels),
    patchFaceRestrictAddressingHost_.setSize(nCreatedLevels),
    meshLevels_.setSize(nCreatedLevels);
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
    const label nCoarseCells
) const
{
    // Check the need for further agglomeration on all processors
    bool contAgg = nCoarseCells >= nCellsInCoarsestLevel_;
    mesh().reduce(contAgg, andOp<bool>());
    return contAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    MeshObject<lduMesh, Foam::GeometricMeshObject, GAMGAgglomeration>(mesh),

    maxLevels_(50),

    nCellsInCoarsestLevel_
    (
        readLabel(controlDict.lookup("nCellsInCoarsestLevel"))
    ),
    meshInterfaces_(mesh.interfaces()),
    nCells_(maxLevels_),
    restrictSortAddressing_(maxLevels_),
    restrictTargetAddressing_(maxLevels_),
    restrictTargetStartAddressing_(maxLevels_),
    restrictAddressingHost_(maxLevels_),
    nFaces_(maxLevels_),
    faceRestrictSortAddressing_(maxLevels_),
    faceRestrictTargetAddressing_(maxLevels_),
    faceRestrictTargetStartAddressing_(maxLevels_),
    faceRestrictAddressingHost_(maxLevels_),
    faceFlipMap_(maxLevels_),
    faceFlipMapHost_(maxLevels_),
    nPatchFaces_(maxLevels_),
    patchFaceRestrictSortAddressing_(maxLevels_),
    patchFaceRestrictTargetAddressing_(maxLevels_),
    patchFaceRestrictTargetStartAddressing_(maxLevels_),
    patchFaceRestrictAddressingHost_(maxLevels_),

    meshLevels_(maxLevels_)
{
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
{
    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            lduMeshConstructorTablePtr_
        );

        lduMeshConstructorTable::iterator cstrIter =
            lduMeshConstructorTablePtr_->find(agglomeratorType);

        if (cstrIter == lduMeshConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "GAMGAgglomeration::New"
                "(const lduMesh& mesh, const dictionary& controlDict)"
            )   << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid matrix GAMGAgglomeration types are :"
                << lduMatrixConstructorTablePtr_->sortedToc() << endl
                << "Valid geometric GAMGAgglomeration types are :"
                << lduMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return store(cstrIter()(mesh, controlDict).ptr());
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
{
    const lduMesh& mesh = matrix.mesh();

    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "algebraicGAMGAgglomerationLibs",
            lduMatrixConstructorTablePtr_
        );

        if
        (
            !lduMatrixConstructorTablePtr_
         || !lduMatrixConstructorTablePtr_->found(agglomeratorType)
        )
        {
            return New(mesh, controlDict);
        }
        else
        {
            lduMatrixConstructorTable::iterator cstrIter =
                lduMatrixConstructorTablePtr_->find(agglomeratorType);

            return store(cstrIter()(matrix, controlDict).ptr());
        }
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


Foam::autoPtr<Foam::GAMGAgglomeration> Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
{
    const word agglomeratorType(controlDict.lookup("agglomerator"));

    const_cast<Time&>(mesh.thisDb().time()).libs().open
    (
        controlDict,
        "geometricGAMGAgglomerationLibs",
        geometryConstructorTablePtr_
    );

    geometryConstructorTable::iterator cstrIter =
        geometryConstructorTablePtr_->find(agglomeratorType);

    if (cstrIter == geometryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "GAMGAgglomeration::New"
            "(const lduMesh& mesh, const scalarField&"
            ", const vectorField&, const dictionary& controlDict)"
        )   << "Unknown GAMGAgglomeration type "
            << agglomeratorType << ".\n"
            << "Valid geometric GAMGAgglomeration types are :"
            << geometryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGAgglomeration>
    (
        cstrIter()
        (
            mesh,
            cellVolumes,
            faceAreas,
            controlDict
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return mesh_;
    }
    else
    {
        return meshLevels_[i - 1];
    }
}


bool Foam::GAMGAgglomeration::hasMeshLevel(const label i) const
{
    if (i == 0)
    {
        return true;
    }
    else
    {
        return meshLevels_.set(i - 1);
    }
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return meshInterfaces_;
    }
    else
    {
        return meshLevels_[i - 1].rawInterfaces();
    }
}


void Foam::GAMGAgglomeration::clearLevel(const label i)
{
    if (hasMeshLevel(i))
    {
        meshLevels_.set(i - 1, NULL);

        if (i < nCells_.size())
        {
            nCells_[i] = -555;
            restrictSortAddressing_.set(i, NULL);
            restrictTargetAddressing_.set(i, NULL);
            restrictTargetStartAddressing_.set(i, NULL);
            restrictAddressingHost_.set(i, NULL);
            nFaces_[i] = -666;
            faceRestrictSortAddressing_.set(i, NULL);
            faceRestrictTargetAddressing_.set(i, NULL);
            faceRestrictTargetStartAddressing_.set(i, NULL);
            faceRestrictAddressingHost_.set(i, NULL);
            faceFlipMap_.set(i, NULL);
            faceFlipMapHost_.set(i, NULL);
            nPatchFaces_.set(i, NULL);
            patchFaceRestrictSortAddressing_.set(i, NULL);
            patchFaceRestrictTargetAddressing_.set(i, NULL);
            patchFaceRestrictTargetStartAddressing_.set(i, NULL);
            patchFaceRestrictAddressingHost_.set(i, NULL);
        }
    }
}


bool Foam::GAMGAgglomeration::checkRestriction
(
    labelList& newRestrict,
    label& nNewCoarse,
    const lduAddressing& fineAddressing,
    const labelUList& restrict,
    const label nCoarse
)
{
    if (fineAddressing.size() != restrict.size())
    {
        FatalErrorIn
        (
            "checkRestriction(..)"
        )   << "nCells:" << fineAddressing.size()
            << " agglom:" << restrict.size()
            << abort(FatalError);
    }

    // Seed (master) for every region
    labelList master(identity(fineAddressing.size()));

    // Now loop and transport master through region
    const labelUList& lower = fineAddressing.lowerAddrHost();
    const labelUList& upper = fineAddressing.upperAddrHost();

    while (true)
    {
        label nChanged = 0;

        forAll(lower, faceI)
        {
            label own = lower[faceI];
            label nei = upper[faceI];

            if (restrict[own] == restrict[nei])
            {
                // coarse-mesh-internal face

                if (master[own] < master[nei])
                {
                    master[nei] = master[own];
                    nChanged++;
                }
                else if (master[own] > master[nei])
                {
                    master[own] = master[nei];
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }


    // Count number of regions/masters per coarse cell
    labelListList coarseToMasters(nCoarse);
    nNewCoarse = 0;
    forAll(restrict, cellI)
    {
        labelList& masters = coarseToMasters[restrict[cellI]];

        if (findIndex(masters, master[cellI]) == -1)
        {
            masters.append(master[cellI]);
            nNewCoarse++;
        }
    }

    if (nNewCoarse > nCoarse)
    {
        //WarningIn("GAMGAgglomeration::checkRestriction(..)")
        //    << "Have " << nCoarse
        //    << " agglomerated cells but " << nNewCoarse
        //    << " disconnected regions" << endl;

        // Keep coarseToMasters[0] the original coarse, allocate new ones
        // for the others
        labelListList coarseToNewCoarse(coarseToMasters.size());

        nNewCoarse = nCoarse;

        forAll(coarseToMasters, coarseI)
        {
            const labelList& masters = coarseToMasters[coarseI];

            labelList& newCoarse = coarseToNewCoarse[coarseI];
            newCoarse.setSize(masters.size());
            newCoarse[0] = coarseI;
            for (label i = 1; i < newCoarse.size(); i++)
            {
                newCoarse[i] = nNewCoarse++;
            }
        }

        newRestrict.setSize(fineAddressing.size());
        forAll(restrict, cellI)
        {
            label coarseI = restrict[cellI];

            label index = findIndex(coarseToMasters[coarseI], master[cellI]);
            newRestrict[cellI] = coarseToNewCoarse[coarseI][index];
        }

        return false;
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
