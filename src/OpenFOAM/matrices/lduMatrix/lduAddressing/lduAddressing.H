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

Class
    Foam::lduAddressing

Description
    The class contains the addressing required by the lduMatrix: upper, lower
    and losort.

    The addressing can be created in two ways: either with references to
    upper and lower in which case it stores references or from labelLists,
    in which case it stores the addressing itself. Additionally, the losort
    addressing belongs to the class is as on lazy evaluation.

    The ordering of owner addresses is such that the labels are in
    increasing order, with groups of identical labels for edges "owned" by
    the same point. The neighbour labels are also ordered in ascending
    order but only for groups of edges belonging to each point. An example
    is given below:
    \verbatim
        owner     neighbour
        0         1
        0         20
        1         2
        1         21
        2         3
        2         22
        3         4
        3         23
        4         5
        4         24
        5         6
        5         25
        6         7
        6         26
        7         8
        7         27
        8         9
        8         28
        9         10
        9         29
    \endverbatim

    There exists an alternative way of addressing the owner
    list: instead of repeating the same label in the owner list, it is
    possible to address the start of each point neighbours in the
    neighbour list. This reduces the size of owner addressing from a list
    over all edges to a list over all points + 1:

    \verbatim
        Owner start list: 0 2 4 6 8 10 12 14 16 18
    \endverbatim

    We shall use the second form of the addressing for fast lookup
    of edge label from the known owner and neighbour, using the following
    algorithm:
    -# take the owner label and position the start of lookup
       using the owner start list
    -# loop through all neighbours of this owner (ending at the start of
      lookup of owner + 1) until the match with current neighbour is found.
      The index used on the neighbour list for the match is the edge index.

    While owner start addressing allows us to find the edge owned by the
    points, it is also necessary to find the edges for which the point is
    a neighbour. Losort addressing lists the edges neighboured by the
    point and we shall use the same trick as above to address into this
    list. Thus, for every point the losort start gives the address of the
    first face to neighbour this point.

SourceFiles
    lduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef lduAddressing_H
#define lduAddressing_H

#include "labelList.H"
#include "lduSchedule.H"
#include "boolList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lduAddressing Declaration
\*---------------------------------------------------------------------------*/

class lduAddressing
{
    // Private data

        //- Number of equations
        label size_;

    //- Demand-driven data

        //- Losort addressing
        mutable labelgpuList* losortPtr_;

        //- Owner start addressing
        mutable labelgpuList* ownerStartPtr_;

        //- Owner sort addressing
        mutable labelgpuList* ownerSortAddrPtr_;

        //- Losort start addressing
        mutable labelgpuList* losortStartPtr_;

        mutable PtrList<const labelgpuList> patchSortCells_;

        mutable PtrList<const labelgpuList> patchSortAddr_;

        mutable PtrList<const labelgpuList> patchSortStartAddr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduAddressing(const lduAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const lduAddressing&);

        //- Calculate losort
        void calcLosort() const;

        //- Calculate owner start
        void calcOwnerStart() const;

        //- Calculate losort start
        void calcLosortStart() const;

        //- Calculate patch sort
        void calcPatchSort() const;

        //- Calculate patch sort start
        void calcPatchSortStart() const;


public:

    // Constructor
    lduAddressing(const label nEqns)
    :
        size_(nEqns),
        losortPtr_(NULL),
        ownerSortAddrPtr_(NULL),
        ownerStartPtr_(NULL),
        losortStartPtr_(NULL)
    {}


    //- Destructor
    virtual ~lduAddressing();


    // Member Functions

        //- Return number of equations
        label size() const
        {
            return size_;
        }

        //- Return number of patches
        virtual label nPatches() const = 0;
        virtual label level() const = 0;

        //- Return lower addressing
        virtual const labelgpuList& lowerAddr() const = 0;
        virtual const labelList& lowerAddrHost() const = 0;

        //- Return upper addressing
        virtual const labelgpuList& upperAddr() const = 0;
        virtual const labelList& upperAddrHost() const = 0;

        //- Return patch to internal addressing given patch number
        virtual bool patchAvailable
        (
            const label patchNo
        ) const = 0;

        virtual const labelgpuList& patchAddr
        (
            const label patchNo
        ) const = 0;

        virtual const labelList& patchAddrHost
        (
            const label patchNo
        ) const = 0;

        const labelgpuList& patchSortCells
        (
            const label patchNo
        ) const;

        const labelgpuList& patchSortAddr
        (
            const label patchNo
        ) const;

        const labelgpuList& patchSortStartAddr
        (
            const label patchNo
        ) const;

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const = 0;

        //- Return losort addressing
        const labelgpuList& losortAddr() const;

        //- Return owner sort addressing
        const labelgpuList& ownerSortAddr() const;

        //- Return owner start addressing
        const labelgpuList& ownerStartAddr() const;

        //- Return losort start addressing
        const labelgpuList& losortStartAddr() const; 

        //- Calculate bandwidth and profile of addressing
        Tuple2<label, scalar> band() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
