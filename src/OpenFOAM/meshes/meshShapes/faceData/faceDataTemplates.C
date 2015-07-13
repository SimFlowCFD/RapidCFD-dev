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

#include "faceData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
__HOST____DEVICE__
Type Foam::faceData::average
(
    const label* labels,
    const point* meshPoints,
    const Type* fld,
    Type zero
) const
{
    // Calculate the average by breaking the face into triangles and
    // area-weighted averaging their averages

    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return
            (1.0/3.0)
           *(
               fld[labels[0]]
             + fld[labels[1]]
             + fld[labels[2]]
            );
    }

    label nPoints = size();

    point centrePoint(0,0,0);
    Type cf = zero;

    for (register label pI=0; pI<nPoints; pI++)
    {
        centrePoint += meshPoints[labels[pI]];
        cf += fld[labels[pI]];
    }

    centrePoint /= nPoints;
    cf /= nPoints;

    scalar sumA = 0;
    Type sumAf = zero;

    for (register label pI=0; pI<nPoints; pI++)
    {
        // Calculate 3*triangle centre field value
        Type ttcf  =
        (
            fld[labels[pI]]
          + fld[labels[(pI + 1) % nPoints]]
          + cf
        );

        // Calculate 2*triangle area
        scalar ta = Foam::mag
        (
            (meshPoints[labels[pI]] - centrePoint)
          ^ (meshPoints[labels[(pI + 1) % nPoints]] - centrePoint)
        );

        sumA += ta;
        sumAf += ta*ttcf;
    }

    if (sumA > VSMALL)
    {
        return sumAf/(3*sumA);
    }
    else
    {
        return cf;
    }
}


// ************************************************************************* //
