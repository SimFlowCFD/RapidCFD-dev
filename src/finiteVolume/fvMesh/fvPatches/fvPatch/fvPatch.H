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
    Foam::fvPatch

Description
    A finiteVolume patch using a polyPatch and a fvBoundaryMesh

SourceFiles
    fvPatch.C
    fvPatchNew.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatch_H
#define fvPatch_H

#include "polyPatch.H"
#include "labelList.H"
#include "SubList.H"
#include "typeInfo.H"
#include "tmp.H"
#include "primitiveFields.H"
#include "SubField.H"
#include "fvPatchFieldsFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvBoundaryMesh;
class surfaceInterpolation;

/*---------------------------------------------------------------------------*\
                           Class fvPatch Declaration
\*---------------------------------------------------------------------------*/

class fvPatch
{
    // Private data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;

        //- Reference to boundary mesh
        const fvBoundaryMesh& boundaryMesh_;


    // Private Member Functions

        //- Disallow construct as copy
        fvPatch(const fvPatch&);

        //- Disallow assignment
        void operator=(const fvPatch&);


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalargpuField&) const;

        //- Initialise the patches for moving points
        virtual void initMovePoints();

        //- Correct patches after moving points
        virtual void movePoints();


public:

    typedef fvBoundaryMesh BoundaryMesh;

    friend class fvBoundaryMesh;
    friend class surfaceInterpolation;

    //- Runtime type information
    TypeName(polyPatch::typeName_());


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvPatch,
            polyPatch,
            (const polyPatch& patch, const fvBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        fvPatch(const polyPatch&, const fvBoundaryMesh&);


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<fvPatch> New
        (
            const polyPatch&,
            const fvBoundaryMesh&
        );


    //- Destructor
    virtual ~fvPatch();


    // Member Functions

    // Access

            //- Return the polyPatch
            const polyPatch& patch() const
            {
                return polyPatch_;
            }

            //- Return name
            const word& name() const
            {
                return polyPatch_.name();
            }

            //- Return start label of this patch in the polyMesh face list
            label start() const
            {
                return polyPatch_.start();
            }

            //- Return size
            virtual label size() const
            {
                return polyPatch_.size();
            }

            //- Return true if this patch is coupled
            virtual bool coupled() const
            {
                return polyPatch_.coupled();
            }

            //- Return true if the given type is a constraint type
            static bool constraintType(const word& pt);

            //- Return a list of all the constraint patch types
            static wordList constraintTypes();

            //- Return the index of this patch in the fvBoundaryMesh
            label index() const
            {
                return polyPatch_.index();
            }

            //- Return boundaryMesh reference
            const fvBoundaryMesh& boundaryMesh() const
            {
                return boundaryMesh_;
            }

            //- Return faceCells
            virtual const labelgpuList& faceCells() const;
            virtual const labelList& faceCellsHost() const;


        // Access functions for geometrical data

            //- Return face centres
            const vectorgpuField& Cf() const;

            //- Return neighbour cell centres
            tmp<vectorgpuField> Cn() const;

            //- Return face area vectors
            const vectorgpuField& Sf() const;

            //- Return face area magnitudes
            const scalargpuField& magSf() const;

            //- Return face normals
            tmp<vectorgpuField> nf() const;

            //- Return cell-centre to face-centre vector
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre vector is returned
            virtual tmp<vectorgpuField> delta() const;


        // Access functions for demand driven data

            //- Return patch weighting factors
            const scalargpuField& weights() const;

            //- Return the face - cell distance coeffient
            //  except for coupled patches for which the cell-centre
            //  to coupled-cell-centre distance coeffient is returned
            const scalargpuField& deltaCoeffs() const;


        // Evaluation functions

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<gpuField<Type> > patchInternalField(const gpuList<Type>&) const;

            //- Return given internal field next to patch as patch field
            template<class Type>
            void patchInternalField(const gpuList<Type>&, gpuField<Type>&) const;

            //- Return given internal field next to patch as patch field
            template<class Type>
            tmp<Field<Type> > patchInternalField(const UList<Type>&) const;

            //- Return given internal field next to patch as patch field
            template<class Type>
            void patchInternalField(const UList<Type>&, Field<Type>&) const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& patchField
            (
                const GeometricField&
            ) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type>
            const typename GeometricField::PatchFieldType& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type* = NULL
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
