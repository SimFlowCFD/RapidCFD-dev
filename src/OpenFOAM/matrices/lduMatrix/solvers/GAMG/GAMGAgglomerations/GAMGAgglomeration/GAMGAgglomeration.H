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
    Foam::GAMGAgglomeration

Description
    Geometric agglomerated algebraic multigrid agglomeration class.

SourceFiles
    GAMGAgglomeration.C
    GAMGAgglomerationTemplates.C
    GAMGAgglomerateLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGAgglomeration_H
#define GAMGAgglomeration_H

#include "MeshObject.H"
#include "lduPrimitiveMesh.H"
#include "lduInterfacePtrsList.H"
#include "primitiveFields.H"
#include "runTimeSelectionTables.H"

#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMesh;
class lduMatrix;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                    Class GAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class GAMGAgglomeration
:
    public MeshObject<lduMesh, GeometricMeshObject, GAMGAgglomeration>
{
protected:

    // Protected data

        //- Max number of levels
        const label maxLevels_;

        //- Number of cells in coarsest level
        const label nCellsInCoarsestLevel_;

        //- Cached mesh interfaces
        const lduInterfacePtrsList meshInterfaces_;

        //- The number of cells in each level
        labelList nCells_;

        //- Cell restriction addressing array.
        //  Maps from the finer to the coarser level.
        PtrList<labelgpuField> restrictSortAddressing_;
        PtrList<labelgpuField> restrictTargetAddressing_;
        PtrList<labelgpuField> restrictTargetStartAddressing_;
        PtrList<labelField> restrictAddressingHost_;

        //- The number of (coarse) faces in each level.
        //  max(faceRestrictAddressing)+1.
        labelList nFaces_;

        //- Face restriction addressing array.
        //  Maps from the finer to the coarser level.
        //  Positive indices map the finer faces which form part of the boundary
        //  of the coarser cells to the corresponding coarser cell face.
        //  Negative indices map the finer faces which are internal to the
        //  coarser cells to minus the corresponding coarser cell index minus 1.
        PtrList<labelgpuList> faceRestrictSortAddressing_;
        PtrList<labelgpuList> faceRestrictTargetAddressing_;
        PtrList<labelgpuList> faceRestrictTargetStartAddressing_;
        PtrList<labelList> faceRestrictAddressingHost_;

        //- Face flip: for faces mapped to internal faces stores whether
        //  the face is reversed or not. This is used to avoid having
        //  to access the coarse mesh at all when mapping
        PtrList<boolgpuList> faceFlipMap_;
        PtrList<boolList> faceFlipMapHost_;

        //- The number of (coarse) patch faces in each level.
        //  max(patchFaceRestrictAddressing_)+1.
        PtrList<labelList> nPatchFaces_;

        //- Patch-local face restriction addressing array.
        //  Maps from the finer to the coarser level. Always positive.
        //  Extracted from GAMGInterfaces after agglomeration.
        PtrList<labelgpuListList> patchFaceRestrictSortAddressing_;
        PtrList<labelgpuListList> patchFaceRestrictTargetAddressing_;
        PtrList<labelgpuListList> patchFaceRestrictTargetStartAddressing_;
        PtrList<labelListList> patchFaceRestrictAddressingHost_;

        //- Hierarchy of mesh addressing
        PtrList<lduPrimitiveMesh> meshLevels_;

    // Protected Member Functions

        //- Assemble coarse mesh addressing
        void agglomerateLduAddressing(const label fineLevelIndex);

        //- Combine a level with the previous one
        void combineLevels(const label curLevel);

        //- Shrink the number of levels to that specified
        void compactLevels(const label nCreatedLevels);

        //- Check the need for further agglomeration
        bool continueAgglomerating(const label nCoarseCells) const;

        void clearLevel(const label leveli);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGAgglomeration(const GAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("GAMGAgglomeration");

       //GPU addressing helper functions
       static void createSort
       (
           const labelgpuList& list, 
           labelgpuList& sort
       );

       static void createTarget
       (
           const labelgpuList& list,
           const labelgpuList& sort,
           labelgpuList& target,
           labelgpuList& targetStart
       );

    // Declare run-time constructor selection tables

        //- Runtime selection table for pure geometric agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGAgglomeration,
            lduMesh,
            (
                const lduMesh& mesh,
                const dictionary& controlDict
            ),
            (
                mesh,
                controlDict
            )
        );

        //- Runtime selection table for matrix or mixed geometric/matrix
        //  agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGAgglomeration,
            lduMatrix,
            (
                const lduMatrix& matrix,
                const dictionary& controlDict
            ),
            (
                matrix,
                controlDict
            )
        );

        //- Runtime selection table for matrix or mixed geometric/matrix
        //  agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGAgglomeration,
            geometry,
            (
                const lduMesh& mesh,
                const scalarField& cellVolumes,
                const vectorField& faceAreas,
                const dictionary& controlDict
            ),
            (
                mesh,
                cellVolumes,
                faceAreas,
                controlDict
            )
        );


    // Constructors

        //- Construct given mesh and controls
        GAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );


    // Selectors

        //- Return the selected geometric agglomerator
        static const GAMGAgglomeration& New
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Return the selected matrix agglomerator
        static const GAMGAgglomeration& New
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );

        //- Return the selected geometric agglomerator
        static autoPtr<GAMGAgglomeration> New
        (
            const lduMesh& mesh,
            const scalarField& cellVolumes,
            const vectorField& faceAreas,
            const dictionary& controlDict
        );


    //- Destructor
    ~GAMGAgglomeration();


    // Member Functions

        // Access

            label size() const
            {
                return meshLevels_.size();
            }

            //- Return LDU mesh of given level
            const lduMesh& meshLevel(const label leveli) const;

            //- Do we have mesh for given level?
            bool hasMeshLevel(const label leveli) const;

            //- Return LDU interface addressing of given level
            const lduInterfacePtrsList& interfaceLevel
            (
                const label leveli
            ) const;

            //- Return cell restrict addressing of given level
            const labelgpuField& restrictSortAddressing(const label leveli) const
            {
                return restrictSortAddressing_[leveli];
            }

            const labelgpuField& restrictTargetAddressing(const label leveli) const
            {
                return restrictTargetAddressing_[leveli];
            }

            const labelgpuField& restrictTargetStartAddressing(const label leveli) const
            {
                return restrictTargetStartAddressing_[leveli];
            }


            const labelField& restrictAddressingHost(const label leveli) const
            {
                return restrictAddressingHost_[leveli];
            }

            //- Face restrict addressing of given level
            const labelgpuList& faceRestrictSortAddressing(const label leveli) const
            {
                return faceRestrictSortAddressing_[leveli];
            }

            const labelgpuList& faceRestrictTargetAddressing(const label leveli) const
            {
                return faceRestrictTargetAddressing_[leveli];
            }

            const labelgpuList& faceRestrictTargetStartAddressing(const label leveli) const
            {
                return faceRestrictTargetStartAddressing_[leveli];
            }


            const labelList& faceRestrictAddressingHost(const label leveli) const
            {
                return faceRestrictAddressingHost_[leveli];
            }


            const labelgpuListList& patchFaceRestrictSortAddressing(const label leveli) const
            {
                return patchFaceRestrictSortAddressing_[leveli];
            }

            const labelgpuListList& patchFaceRestrictTargetAddressing(const label leveli) const
            {
                return patchFaceRestrictTargetAddressing_[leveli];
            }

            const labelgpuListList& patchFaceRestrictTargetStartAddressing(const label leveli) const
            {
                return patchFaceRestrictTargetStartAddressing_[leveli];
            }


            const labelListList& patchFaceRestrictAddressingHost(const label leveli) const
            {
                return patchFaceRestrictAddressingHost_[leveli];
            }

            //- Return face flip map of given level
            const boolgpuList& faceFlipMap(const label leveli) const
            {
                return faceFlipMap_[leveli];
            }

            const boolList& faceFlipMapHost(const label leveli) const
            {
                return faceFlipMapHost_[leveli];
            }

            //- Return number of coarse cells
            label nCells(const label leveli) const
            {
                return nCells_[leveli];
            }

            //- Return number of coarse faces
            label nFaces(const label leveli) const
            {
                return nFaces_[leveli];
            }

            //- Return number of coarse patch faces
            const labelList& nPatchFaces(const label leveli) const
            {
                return nPatchFaces_[leveli];
            }


        // Restriction and prolongation

            //- Restrict (integrate by summation) cell field
            template<class Type>
            void restrictField
            (
                gpuField<Type>& cf,
                const gpuField<Type>& ff,
                const label fineLevelIndex
            ) const;

            //- Restrict (integrate by summation) face field
            template<class Type>
            void restrictFaceField
            (
                gpuField<Type>& cf,
                const gpuField<Type>& ff,
                const label fineLevelIndex
            ) const;

            template<class Type>
            void restrictFaceField
            (
                Field<Type>& cf,
                const Field<Type>& ff,
                const label fineLevelIndex
            ) const;

            //- Restrict (integrate by summation) cell field
            template<class Type>
            void restrictField
            (
                gpuField<Type>& cf,
                const gpuField<Type>& ff,
                const labelgpuList& sort,
                const labelgpuList& target,
                const labelgpuList& targetStart
            ) const;


            //- Prolong (interpolate by injection) cell field
            template<class Type>
            void prolongField
            (
                gpuField<Type>& ff,
                const gpuField<Type>& cf,
                const label coarseLevelIndex
            ) const;

        //- Given restriction determines if coarse cells are connected.
        //  Return ok is so, otherwise creates new restriction that is
        static bool checkRestriction
        (
            labelList& newRestrict,
            label& nNewCoarse,
            const lduAddressing& fineAddressing,
            const labelUList& restrict,
            const label nCoarse
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GAMGAgglomerationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
