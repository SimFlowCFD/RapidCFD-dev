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
    Foam::GAMGInterface

Description
     Abstract base class for GAMG agglomerated interfaces.

SourceFiles
    GAMGInterface.C
    newAmgInterface.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGInterface_H
#define GAMGInterface_H

#include "autoPtr.H"
#include "lduInterfacePtrsList.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class GAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class GAMGInterface
:
    public lduInterface
{

protected:

    // Protected data

        //- My index in coarseInterfaces
        const label index_;

        //- All interfaces
        const lduInterfacePtrsList& coarseInterfaces_;

        //- Face-cell addressing
        labelgpuList faceCells_;
        labelgpuList sortCells_;
        labelgpuList cellFaces_;
        labelgpuList cellFacesStart_;
        labelList faceCellsHost_;

        //- Face restrict addressing
        labelgpuList faceRestrictSortAddressing_;
        labelgpuList faceRestrictTargetAddressing_;
        labelgpuList faceRestrictTargetStartAddressing_;
        labelList faceRestrictAddressingHost_;

        void updateAddressing();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGInterface(const GAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGInterface&);


public:

    //- Runtime type information
    TypeName("GAMGInterface");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterface,
            lduInterface,
            (
                const label index,
                const lduInterfacePtrsList& coarseInterfaces,
                const lduInterface& fineInterface,
                const labelField& localRestrictAddressing,
                const labelField& neighbourRestrictAddressing,
                const label fineLevelIndex,
                const label coarseComm
            ),
            (
                index,
                coarseInterfaces,
                fineInterface,
                localRestrictAddressing,
                neighbourRestrictAddressing,
                fineLevelIndex,
                coarseComm
            )
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterface,
            Istream,
            (
                const label index,
                const lduInterfacePtrsList& coarseInterfaces,
                Istream& is
            ),
            (
                index,
                coarseInterfaces,
                is
            )
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterface> New
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterface> New
        (
            const word& coupleType,
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    // Constructors

        //- Construct from interfaces, restrict addressing set later on
        GAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces
        )
        :
            index_(index),
            coarseInterfaces_(coarseInterfaces)
        {}


        //- Construct from interfaces and restrict addressing
        GAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddressing
        )
        :
            index_(index),
            coarseInterfaces_(coarseInterfaces),
            faceCellsHost_(faceCells),
            faceCells_(faceCellsHost_),
            faceRestrictAddressingHost_(faceRestrictAddressing)
        {
            updateAddressing();
        }


        //- Construct from Istream
        GAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    // Member Functions

        // Access

            //- Return size
            virtual label size() const
            {
                return faceCells_.size();
            }

            virtual label index() const
            {
                return index_;
            }

            virtual const lduInterfacePtrsList& coarseInterfaces() const
            {
                return coarseInterfaces_;
            }

            //- Return faceCell addressing
            virtual const labelgpuList& faceCells() const
            {
                return faceCells_;
            }

            virtual const labelgpuList& sortCells() const
            {
                return sortCells_;
            }

            virtual const labelgpuList& cellFaces() const
            {
                return cellFaces_;
            }

            virtual const labelgpuList& cellFacesStart() const
            {
                return cellFacesStart_;
            }

            virtual const labelList& faceCellsHost() const
            {
                return faceCellsHost_;
            }

            //- Return (local)face restrict addressing
            virtual const labelgpuList& faceRestrictSortAddressing() const
            {
                return faceRestrictSortAddressing_;
            }

            virtual const labelgpuList& faceRestrictTargetAddressing() const
            {
                return faceRestrictTargetAddressing_;
            }

            virtual const labelgpuList& faceRestrictTargetStartAddressing() const
            {
                return faceRestrictTargetStartAddressing_;
            }

            virtual const labelList& faceRestrictAddressingHost() const
            {
                return faceRestrictAddressingHost_;
            }

            //- Return the interface internal field of the given field
            template<class Type>
            tmp<gpuField<Type> > interfaceInternalField
            (
                const gpuList<Type>& internalData
            ) const;

            //- Get the interface internal field of the given field
            template<class Type>
            void interfaceInternalField
            (
                const gpuList<Type>& internalData,
                gpuList<Type>&
            ) const;

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;


        // Agglomeration

            //- Merge the next level with this level
            //  combining the face-restrict addressing
            //  and copying the face-cell addressing
            void combine(const GAMGInterface&);

            //- Agglomerating the given fine-level coefficients and return
            virtual tmp<scalargpuField> agglomerateCoeffs
            (
                const scalargpuField& fineCoeffs
            ) const;


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
