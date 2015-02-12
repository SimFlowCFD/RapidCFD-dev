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

Class
    Foam::cellData

Description
    A cellD is defined as a list of faces with extra functionality

SourceFiles
    cellDataI.H

\*---------------------------------------------------------------------------*/

#ifndef cellData_H
#define cellData_H

#include "faceList.H"
#include "oppositeFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cellData;
__HOST____DEVICE__
bool operator==(const cellData&, const cellData&);
__HOST____DEVICE__
inline bool operator!=(const cellData&, const cellData&);

typedef gpuList<cellData> cellDatagpuList;
typedef List<cellData> cellDataList;


/*---------------------------------------------------------------------------*\
                           Class cellData Declaration
\*---------------------------------------------------------------------------*/

class cellData
{
private:
    label start;
    label size;

public:

    // Static data members

        static const char* const typeName;
        static const cellData zero;


    // Constructors

        //- Construct null
        __HOST____DEVICE__
        inline cellData();

        //- Construct given size
        __HOST____DEVICE__
        explicit inline cellData(label,label);

    // Member Functions

        //- Return number of faces
        __HOST____DEVICE__
        inline label nFaces() const
        {
            return size;
        }

        __HOST____DEVICE__
        inline label getStart() const
        {
            return start;
        }

        //- Return labels of cell vertices
//        labelList labels(const faceUList&) const;

        //- Return the cell vertices
//        pointField points(const faceUList&, const pointField&) const;

        //- Return cell edges
//        edgeList edges(const faceUList&) const;

        //- Return index of opposite face
/*        label opposingFaceLabel
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;

        //- Return opposite face oriented the same way as the master face
        oppositeFace opposingFace
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;
*/

        // Instructions for centre and mag:
        // When one wants to access the cell centre and magnitude, the
        // functionality on the mesh level should be used in preference to the
        // functions provided here. They do not rely to the functionality
        // implemented here, provide additional checking and are more efficient.
        // The cell::centre and cell::mag functions may be removed in the
        // future.

        //- Returns cell centre
//        __HOST____DEVICE__
//        point centre(const pointField&, const faceUList&) const;

        //- Returns cell volume
//        __HOST____DEVICE__
//        scalar mag(const pointField&, const faceUList&) const;


    // Friend Operators

        friend bool operator==(const cellData&, const cellData&);
        friend bool operator!=(const cellData&, const cellData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
