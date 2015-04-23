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
    Foam::GAMGSolver

Description
    Geometric agglomerated algebraic multigrid solver.

  Characteristics:
      - Requires positive definite, diagonally dominant matrix.
      - Agglomeration algorithm: selectable and optionally cached.
      - Restriction operator: summation.
      - Prolongation operator: injection.
      - Smoother: Gauss-Seidel.
      - Coarse matrix creation: central coefficient: summation of fine grid
        central coefficients with the removal of intra-cluster face;
        off-diagonal coefficient: summation of off-diagonal faces.
      - Coarse matrix scaling: performed by correction scaling, using steepest
        descent optimisation.
      - Type of cycle: V-cycle with optional pre-smoothing.
      - Coarsest-level matrix solved using ICCG or BICCG.

SourceFiles
    GAMGSolver.C
    GAMGSolverAgglomerateMatrix.C
    GAMGSolverInterpolate.C
    GAMGSolverScale.C
    GAMGSolverSolve.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGSolver_H
#define GAMGSolver_H

#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GAMGSolver Declaration
\*---------------------------------------------------------------------------*/

class GAMGSolver
:
    public lduMatrix::solver
{
    // Private data

        bool cacheAgglomeration_;

        //- Number of pre-smoothing sweeps
        label nPreSweeps_;

        //- Lever multiplier for the number of pre-smoothing sweeps
        label preSweepsLevelMultiplier_;

        //- Maximum number of pre-smoothing sweeps
        label maxPreSweeps_;

        //- Number of post-smoothing sweeps
        label nPostSweeps_;

        //- Lever multiplier for the number of post-smoothing sweeps
        label postSweepsLevelMultiplier_;

        //- Maximum number of post-smoothing sweeps
        label maxPostSweeps_;

        //- Number of smoothing sweeps on finest mesh
        label nFinestSweeps_;

        //- Choose if the corrections should be interpolated after injection.
        //  By default corrections are not interpolated.
        bool interpolateCorrection_;

        //- Choose if the corrections should be scaled.
        //  By default corrections for symmetric matrices are scaled
        //  but not for asymmetric matrices.
        bool scaleCorrection_;

        //- The agglomeration
        const GAMGAgglomeration& agglomeration_;

        //- Hierarchy of matrix levels
        PtrList<lduMatrix> matrixLevels_;

        //- Hierarchy of interfaces.
        PtrList<PtrList<lduInterfaceField> > primitiveInterfaceLevels_;

        //- Hierarchy of interfaces in lduInterfaceFieldPtrs form
        PtrList<lduInterfaceFieldPtrsList> interfaceLevels_;

        //- Hierarchy of interface boundary coefficients
        PtrList<FieldField<gpuField, scalar> > interfaceLevelsBouCoeffs_;

        //- Hierarchy of interface internal coefficients
        PtrList<FieldField<gpuField, scalar> > interfaceLevelsIntCoeffs_;


    // Private Member Functions

        //- Read control parameters from the control dictionary
        virtual void readControls();

        //- Simplified access to interface level
        const lduInterfaceFieldPtrsList& interfaceLevel
        (
            const label i
        ) const;

        //- Simplified access to matrix level
        const lduMatrix& matrixLevel(const label i) const;

        //- Simplified access to interface boundary coeffs level
        const FieldField<gpuField, scalar>& interfaceBouCoeffsLevel
        (
            const label i
        ) const;

        //- Simplified access to interface internal coeffs level
        const FieldField<gpuField, scalar>& interfaceIntCoeffsLevel
        (
            const label i
        ) const;

        //- Agglomerate coarse matrix. Supply mesh to use - so we can
        //  construct temporary matrix on the fine mesh (instead of the coarse
        //  mesh)
        void agglomerateMatrix
        (
            const label fineLevelIndex,
            const lduMesh& coarseMesh,
            const lduInterfacePtrsList& coarseMeshInterfaces
        );

        //- Agglomerate coarse interface coefficients
        void agglomerateInterfaceCoefficients
        (
            const label fineLevelIndex,
            const lduInterfacePtrsList& coarseMeshInterfaces,
            PtrList<lduInterfaceField>& coarsePrimInterfaces,
            lduInterfaceFieldPtrsList& coarseInterfaces,
            FieldField<gpuField, scalar>& coarseInterfaceBouCoeffs,
            FieldField<gpuField, scalar>& coarseInterfaceIntCoeffs
        ) const;

        //- Interpolate the correction after injected prolongation
        void interpolate
        (
            scalargpuField& psi,
            scalargpuField& Apsi,
            const lduMatrix& m,
            const FieldField<gpuField, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt
        ) const;

        //- Interpolate the correction after injected prolongation and
        /// re-normalise
        void interpolate
        (
            scalargpuField& psi,
            scalargpuField& Apsi,
            const lduMatrix& m,
            const FieldField<gpuField, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const labelgpuList& restrictSortAddressing,
            const labelgpuList& restrictTargetAddressing,
            const labelgpuList& restrictTargetStartAddressing,
            const scalargpuField& psiC,
            const direction cmpt
        ) const;

        //- Calculate and apply the scaling factor from Acf, coarseSource
        //  and coarseField.
        //  At the same time do a Jacobi iteration on the coarseField using
        //  the Acf provided after the coarseField values are used for the
        //  scaling factor.
        void scale
        (
            scalargpuField& field,
            scalargpuField& Acf,
            const lduMatrix& A,
            const FieldField<gpuField, scalar>& interfaceLevelBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaceLevel,
            const scalargpuField& source,
            const direction cmpt
        ) const;

        //- Initialise the data structures for the V-cycle
        void initVcycle
        (
            PtrList<scalargpuField>& coarseCorrFields,
            PtrList<scalargpuField>& coarseSources,
            PtrList<lduMatrix::smoother>& smoothers,
            scalargpuField& scratch1,
            scalargpuField& scratch2
        ) const;


        //- Perform a single GAMG V-cycle with pre, post and finest smoothing.
        void Vcycle
        (
            const PtrList<lduMatrix::smoother>& smoothers,
            scalargpuField& psi,
            const scalargpuField& source,
            scalargpuField& Apsi,
            scalargpuField& finestCorrection,
            scalargpuField& finestResidual,

            scalargpuField& scratch1,
            scalargpuField& scratch2,

            PtrList<scalargpuField>& coarseCorrFields,
            PtrList<scalargpuField>& coarseSources,
            const direction cmpt=0
        ) const;


        //- Solve the coarsest level with either an iterative or direct solver
        void solveCoarsestLevel
        (
            scalargpuField& coarsestCorrField,
            const scalargpuField& coarsestSource
        ) const;


public:

    friend class GAMGPreconditioner;

    //- Runtime type information
    TypeName("GAMG");


    // Constructors

        //- Construct from lduMatrix and solver controls
        GAMGSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<gpuField, scalar>& interfaceBouCoeffs,
            const FieldField<gpuField, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~GAMGSolver();


    // Member Functions

        //- Solve
        virtual solverPerformance solve
        (
            scalargpuField& psi,
            const scalargpuField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
