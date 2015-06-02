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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

__host__ __device__
Foam::tetIndices::tetIndices()
:
    cellI_(-1),
    faceI_(-1),
    faceBasePtI_(-1),
    facePtAI_(-1),
    facePtBI_(-1),
    tetPtI_(-1)
{}

__host__ __device__
Foam::tetIndices::tetIndices
(
    label cellI,
    label faceI,
    label faceBasePtI,
    label facePtAI,
    label facePtBI,
    label tetPtI
)
:
    cellI_(cellI),
    faceI_(faceI),
    faceBasePtI_(faceBasePtI),
    facePtAI_(facePtAI),
    facePtBI_(facePtBI),
    tetPtI_(tetPtI)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

__host__ __device__
Foam::label Foam::tetIndices::cell() const
{
    return cellI_;
}


__host__ __device__
Foam::label Foam::tetIndices::face() const
{
    return faceI_;
}


__host__ __device__
Foam::label Foam::tetIndices::faceBasePt() const
{
    return faceBasePtI_;
}


__host__ __device__
Foam::label Foam::tetIndices::facePtA() const
{
    return facePtAI_;
}


__host__ __device__
Foam::label Foam::tetIndices::facePtB() const
{
    return facePtBI_;
}


__host__ __device__
Foam::label Foam::tetIndices::tetPt() const
{
    return tetPtI_;
}


Foam::tetPointRef Foam::tetIndices::tet(const polyMesh& mesh) const
{
    const pointField& pPts = mesh.points();
    const faceList& pFaces = mesh.faces();
    const vectorField& pC = mesh.cellCentres();

    const Foam::face& f = pFaces[faceI_];

    return tetPointRef
    (
        pC[cellI_],
        pPts[f[faceBasePtI_]],
        pPts[f[facePtAI_]],
        pPts[f[facePtBI_]]
    );
}


Foam::tetPointRef Foam::tetIndices::oldTet(const polyMesh& mesh) const
{
    const pointField& oldPPts = mesh.oldPoints();
    const faceList& pFaces = mesh.faces();

    // We need to reconstruct the old Cc from oldPoints (it isn't
    // stored)
    point oldC = mesh.cells()[cellI_].centre
    (
        oldPPts,
        pFaces
    );

    const Foam::face& f = pFaces[faceI_];

    return tetPointRef
    (
        oldC,
        oldPPts[f[faceBasePtI_]],
        oldPPts[f[facePtAI_]],
        oldPPts[f[facePtBI_]]
    );
}


Foam::triPointRef Foam::tetIndices::faceTri(const polyMesh& mesh) const
{
    const pointField& pPts = mesh.points();
    const faceList& pFaces = mesh.faces();

    const Foam::face& f = pFaces[faceI_];

    return triPointRef
    (
        pPts[f[faceBasePtI_]],
        pPts[f[facePtAI_]],
        pPts[f[facePtBI_]]
    );
}


Foam::triFace Foam::tetIndices::faceTriIs(const polyMesh& mesh) const
{
    const faceList& pFaces = mesh.faces();

    const Foam::face& f = pFaces[faceI_];

    return triFace
    (
        f[faceBasePtI_],
        f[facePtAI_],
        f[facePtBI_]
    );
}


Foam::triPointRef Foam::tetIndices::oldFaceTri(const polyMesh& mesh) const
{
    const pointField& oldPPts = mesh.oldPoints();
    const faceList& pFaces = mesh.faces();

    const Foam::face& f = pFaces[faceI_];

    return triPointRef
    (
        oldPPts[f[faceBasePtI_]],
        oldPPts[f[facePtAI_]],
        oldPPts[f[facePtBI_]]
    );
}


__host__ __device__
Foam::label& Foam::tetIndices::cell()
{
    return cellI_;
}


__host__ __device__
Foam::label& Foam::tetIndices::face()
{
    return faceI_;
}


__host__ __device__
Foam::label& Foam::tetIndices::faceBasePt()
{
    return faceBasePtI_;
}


__host__ __device__
Foam::label& Foam::tetIndices::facePtA()
{
    return facePtAI_;
}


__host__ __device__
Foam::label& Foam::tetIndices::facePtB()
{
    return facePtBI_;
}


__host__ __device__
Foam::label& Foam::tetIndices::tetPt()
{
    return tetPtI_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

__host__ __device__
inline bool Foam::tetIndices::operator==(const Foam::tetIndices& rhs) const
{
    return
    (
        cell() == rhs.cell()
     && face() == rhs.face()
     && faceBasePt() == rhs.faceBasePt()
     && facePtA() == rhs.facePtA()
     && facePtB() == rhs.facePtB()
     && tetPt() == rhs.tetPt()
    );
}


__host__ __device__
inline bool Foam::tetIndices::operator!=(const Foam::tetIndices& rhs) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
