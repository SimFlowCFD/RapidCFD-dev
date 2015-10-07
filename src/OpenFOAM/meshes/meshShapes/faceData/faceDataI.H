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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
__HOST____DEVICE__
inline Foam::faceData::faceData()
:
   start_(0),
   size_(0)
{}

__HOST____DEVICE__
inline Foam::faceData::faceData(label start, label size)
:
   start_(start),
   size_(size)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline Foam::scalar Foam::faceData::mag(const label* labels,const point* p) const
{
    return ::Foam::mag(normal(labels,p));
}

// Next vertex on faceData
__HOST____DEVICE__
inline Foam::label Foam::faceData::nextLabel(const label* labels,const label i) const
{
    return labels[start()+i+1];
}


// Previous vertex on faceData
__HOST____DEVICE__
inline Foam::label Foam::faceData::prevLabel(const label* labels,const label i) const
{
    return labels[start()+i-1];
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
__HOST____DEVICE__
inline bool Foam::operator==(const faceData& a, const faceData& b)
{
    return a.size() == b.size() && a.start() == b.start();
}

__HOST____DEVICE__
inline bool Foam::operator!=(const faceData& a, const faceData& b)
{
    return ! operator==(a,b);
}

__HOST____DEVICE__
inline Foam::point Foam::faceData::centre(const label* labels, const point* points) const
{
    // Calculate the centre by breaking the faceData into triangles and
    // area-weighted averaging their centres

    const label nPoints = size();

    // If the faceData is a triangle, do a direct calculation
    if (nPoints == 3)
    {
        return
            (1.0/3.0)
           *(
               points[labels[0]]
             + points[labels[1]]
             + points[labels[2]]
            );
    }


    point centrePoint(0,0,0);
    for (register label pI=0; pI<nPoints; ++pI)
    {
        centrePoint += points[labels[pI]];
    }
    centrePoint /= nPoints;

    scalar sumA = 0;
    vector sumAc(0,0,0);

    for (register label pI=0; pI<nPoints; ++pI)
    {
        const point& nextPoint = points[labels[(pI + 1) % nPoints]];

        // Calculate 3*triangle centre
        const vector ttc
        (
            points[labels[pI]]
          + nextPoint
          + centrePoint
        );

        // Calculate 2*triangle area
        const scalar ta = Foam::mag
        (
            (points[labels[pI]] - centrePoint)
          ^ (nextPoint - centrePoint)
        );

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        return sumAc/(3.0*sumA);
    }
    else
    {
        return centrePoint;
    }
}

__HOST____DEVICE__
inline Foam::vector Foam::faceData::normal(const label* labels, const point* p) const
{
    const label nPoints = size();

    // Calculate the normal by summing the faceData triangle normals.
    // Changed to deal with small concavity by using a central decomposition
    //

    // If the faceData is a triangle, do a direct calculation to avoid round-off
    // error-related problems
    //
    if (nPoints == 3)
    {
        return triPointRef
        (
            p[labels[0]],
            p[labels[1]],
            p[labels[2]]
        ).normal();
    }

    register label pI;

    point centrePoint(0,0,0);
    for (pI = 0; pI < nPoints; ++pI)
    {
        centrePoint += p[labels[pI]];
    }
    centrePoint /= nPoints;

    vector n(0,0,0);

    point nextPoint = centrePoint;

    for (pI = 0; pI < nPoints; ++pI)
    {
        if (pI < nPoints - 1)
        {
            nextPoint = p[labels[pI + 1]];
        }
        else
        {
            nextPoint = p[labels[0]];
        }

        // Note: for best accuracy, centre point always comes last
        //
        n += triPointRef
        (
            p[labels[pI]],
            nextPoint,
            centrePoint
        ).normal();
    }

    return n;
}

__HOST____DEVICE__
inline Foam::scalar Foam::faceData::sweptVol
(
    const label* labels,
    const point* oldPoints,
    const point* newPoints
) const
{
    // This Optimization causes a small discrepancy between the swept-volume of
    // opposite faceDatas of complex cells with triangular faceDatas opposing polygons.
    // It could be used without problem for tetrahedral cells
    // if (size() == 3)
    // {
    //     return
    //     (
    //         triPointRef
    //         (
    //             oldPoints[operator[](0)],
    //             oldPoints[operator[](1)],
    //             oldPoints[operator[](2)]
    //         ).sweptVol
    //         (
    //             triPointRef
    //             (
    //                 newPoints[operator[](0)],
    //                 newPoints[operator[](1)],
    //                 newPoints[operator[](2)]
    //             )
    //         )
    //     );
    // }

    scalar sv = 0;

    // Calculate the swept volume by breaking the faceData into triangles and
    // summing their swept volumes.
    // Changed to deal with small concavity by using a central decomposition

    point centreOldPoint = centre(labels,oldPoints);
    point centreNewPoint = centre(labels,newPoints);

    label nPoints = size();

    for (register label pi=0; pi<nPoints-1; ++pi)
    {
        // Note: for best accuracy, centre point always comes last
        sv += triPointRef
        (
            centreOldPoint,
            oldPoints[labels[pi]],
            oldPoints[labels[pi + 1]]
        ).sweptVol
        (
            triPointRef
            (
                centreNewPoint,
                newPoints[labels[pi]],
                newPoints[labels[pi + 1]]
            )
        );
    }

    sv += triPointRef
    (
        centreOldPoint,
        oldPoints[labels[nPoints-1]],
        oldPoints[labels[0]]
    ).sweptVol
    (
        triPointRef
        (
            centreNewPoint,
            newPoints[labels[nPoints-1]],
            newPoints[labels[0]]
        )
    );

    return sv;
}

__HOST____DEVICE__
inline Foam::scalar Foam::faceData::contactSphereDiameter
(
    const label* labels,
    const point& p,
    const vector& n,
    const point* meshPoints
) const
{
    scalar magN = Foam::mag(n);

    vector n1 = n/(magN + SMALL);
    vector n2 = normal(labels,meshPoints);

    n2 /= Foam::mag(n2) + SMALL;

    return 2*((centre(labels,meshPoints) - p) & n2)/((n1 & n2) - 1.0);
}


// ************************************************************************* //
