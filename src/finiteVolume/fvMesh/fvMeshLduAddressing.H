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
    Foam::fvMeshLduAddressing

Description
    Foam::fvMeshLduAddressing

SourceFiles
    fvMeshLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshLduAddressing_H
#define fvMeshLduAddressing_H

#include "lduAddressing.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class fvMeshLduAddressing
:
    public lduAddressing
{
    // Private data

        //- Lower as a subList of allOwner
        labelgpuList lowerAddr_;

        //- Upper as a reference to neighbour
        const labelgpuList& upperAddr_;

        //- Patch addressing as a list of sublists
        List<const labelgpuList*> patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;

        //- Lower as a subList of allOwner
        labelList::subList lowerAddrHost_;

        //- Upper as a reference to neighbour
        const labelList& upperAddrHost_;

        //- Patch addressing as a list of sublists
        List<const labelList*> patchAddrHost_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMeshLduAddressing(const fvMeshLduAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshLduAddressing&);


public:

    // Constructors

        //- Construct from components
        fvMeshLduAddressing(const fvMesh& mesh)
        :
            lduAddressing(mesh.nCells()),
            lowerAddr_
            (
                mesh.getFaceOwner(),
                mesh.nInternalFaces()
            ),
            upperAddr_(mesh.getFaceNeighbour()),
            patchAddr_(mesh.boundary().size()),
            lowerAddrHost_
            (
                labelList::subList
                (
                    mesh.faceOwner(),
                    mesh.nInternalFaces()
                )
            ),
            upperAddrHost_(mesh.faceNeighbour()),
            patchAddrHost_(mesh.boundary().size()),
            patchSchedule_(mesh.globalData().patchSchedule())
        {
            forAll(mesh.boundary(), patchI)
            {
                patchAddr_[patchI] = &mesh.boundary()[patchI].faceCells();
            }
            forAll(mesh.boundary(), patchI)
            {
                patchAddrHost_[patchI] = &mesh.boundary()[patchI].faceCellsHost();
            }
        }


    //- Destructor
    ~fvMeshLduAddressing()
    {}


    // Member Functions
        
        
        label nPatches() const
        {
            return patchAddr_.size();
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        const labelgpuList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        const labelgpuList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return wheather patch addressing is available
        bool patchAvailable(const label i) const
        {
            return true;
        }

        label level() const
        {
            return 0;
        }

        const labelgpuList& patchAddr(const label i) const
        {
            return *patchAddr_[i];
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        const labelList& lowerAddrHost() const
        {
            return lowerAddrHost_;
        }

        //- Return upper addressing (i.e. upper label)
        const labelList& upperAddrHost() const
        {
            return upperAddrHost_;
        }

        //- Return patch addressing
        const labelList& patchAddrHost(const label i) const
        {
            return *patchAddrHost_[i];
        }

        // Return patch field evaluation schedule
        const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
