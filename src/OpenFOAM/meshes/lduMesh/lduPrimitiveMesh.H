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
    Foam::lduPrimitiveMesh

Description
    Simplest contrete lduMesh which stores the addressing needed by lduMatrix.

SourceFiles
    lduPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMesh_H
#define lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private data
        label level_;

        //- Lower addressing
        labelgpuList lowerAddr_;
        labelList lowerAddrHost_;

        //- Upper addressing
        labelgpuList upperAddr_;
        labelList upperAddrHost_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces_;

        //- Concrete interfaces
        PtrList<const lduInterface> primitiveInterfaces_;

        //- Patch field evaluation schedule
        lduSchedule patchSchedule_;

        //- Communicator to use for any parallel communication
        const label comm_;

    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const PtrList<lduPrimitiveMesh>&);

        static labelList upperTriOrder
        (
            const label nCells,
            const labelUList& lower,
            const labelUList& upper
        );

        //- Check if in upper-triangular ordering
        static void checkUpperTriangular
        (
            const label size,
            const labelUList& l,
            const labelUList& u
        );

        //- Disallow default bitwise copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMesh&);


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("lduPrimitiveMesh");

    // Constructors

        //- Construct from components but without interfaces. Add interfaces
        //  separately using addInterfaces
        lduPrimitiveMesh
        (
            const label level,
            const label nCells,
            labelgpuList& l,
            labelgpuList& u,
            const label comm,
            bool reUse
        );

        lduPrimitiveMesh
        (
            const label level,
            const label nCells,
            labelList& l,
            labelList& u,
            const label comm,
            bool reUse
        );

        //- Add interfaces to a mesh constructed without
        void addInterfaces
        (
            lduInterfacePtrsList& interfaces,
            const lduSchedule& ps
        );

        //- Construct from components and re-use storage.
        lduPrimitiveMesh
        (
            const label level,
            const label nCells,
            labelList& l,
            labelList& u,
            PtrList<const lduInterface>& primitiveInterfaces,
            const lduSchedule& ps,
            const label comm
        );

        //- Construct from components and re-use storage.
        lduPrimitiveMesh
        (
            const label level,
            const label nCells,
            labelgpuList& l,
            labelgpuList& u,
            PtrList<const lduInterface>& primitiveInterfaces,
            const lduSchedule& ps,
            const label comm
        );

    //- Destructor
    virtual ~lduPrimitiveMesh()
    {}


    // Member Functions

        // Access

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            virtual label nPatches() const
            {
                return interfaces_.size();
            }

            virtual label level() const
            {
                return level_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            //  (reference to cached interfaces)
            const lduInterfacePtrsList& rawInterfaces() const
            {
                return interfaces_;
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return comm_;
            }

            //- Return Lower addressing
            virtual const labelgpuList& lowerAddr() const
            {
                return lowerAddr_;
            }

            virtual const labelList& lowerAddrHost() const
            {
                return lowerAddrHost_;
            }

            //- Return Upper addressing
            virtual const labelgpuList& upperAddr() const
            {
                return upperAddr_;
            }

            virtual const labelList& upperAddrHost() const
            {
                return upperAddrHost_;
            }

            //- Return patch addressing
            virtual bool patchAvailable(const label i) const
            {
                return interfaces_.set(i);
            }

            virtual const labelgpuList& patchAddr(const label i) const
            {
                return interfaces_[i].faceCells();
            }

            virtual const labelList& patchAddrHost(const label i) const
            {
                notImplemented("patchAddrHost(const label i) const");
                return labelList::null();
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }


        // Helper

            //- Get non-scheduled send/receive schedule
            template<class ProcPatch>
            static lduSchedule nonBlockingSchedule(const lduInterfacePtrsList&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lduPrimitiveMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
