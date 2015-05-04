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

\*---------------------------------------------------------------------------*/

#include "GAMGSolver.H"
#include "GAMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<GAMGSolver>
        addGAMGSolverMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<GAMGSolver>
        addGAMGAsymSolverMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGSolver::GAMGSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const FieldField<gpuField, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),

    cacheAgglomeration_(true),
    nPreSweeps_(0),
    preSweepsLevelMultiplier_(1),
    maxPreSweeps_(4),
    nPostSweeps_(2),
    postSweepsLevelMultiplier_(1),
    maxPostSweeps_(4),
    nFinestSweeps_(2),
    interpolateCorrection_(false),
    scaleCorrection_(matrix.symmetric()),
    agglomeration_(GAMGAgglomeration::New(matrix_, controlDict_)),

    matrixLevels_(agglomeration_.size()),
    primitiveInterfaceLevels_(agglomeration_.size()),
    interfaceLevels_(agglomeration_.size()),
    interfaceLevelsBouCoeffs_(agglomeration_.size()),
    interfaceLevelsIntCoeffs_(agglomeration_.size())
{
    readControls();

    forAll(agglomeration_, fineLevelIndex)
    {
        // Agglomerate on to coarse level mesh
        agglomerateMatrix
        (
            fineLevelIndex,
            agglomeration_.meshLevel(fineLevelIndex + 1),
            agglomeration_.interfaceLevel(fineLevelIndex + 1)
        );
    }

    // Tell coarsest matrix its status
    if(agglomeration_.size())
    {
        matrixLevels_[agglomeration_.size()-1].coarsestLevel() = true;
    }

    if (debug)
    {
        for
        (
            label fineLevelIndex = 0;
            fineLevelIndex <= matrixLevels_.size();
            fineLevelIndex++
        )
        {
            if (fineLevelIndex == 0 || matrixLevels_.set(fineLevelIndex-1))
            {
                const lduMatrix& matrix = matrixLevel(fineLevelIndex);
                const lduInterfaceFieldPtrsList& interfaces =
                    interfaceLevel(fineLevelIndex);

                Pout<< "level:" << fineLevelIndex << nl
                    << "    nCells:" << matrix.diag().size() << nl
                    << "    nFaces:" << matrix.lower().size() << nl
                    << "    nInterfaces:" << interfaces.size()
                    << endl;

                forAll(interfaces, i)
                {
                    if (interfaces.set(i))
                    {
                        Pout<< "        " << i
                            << "\ttype:" << interfaces[i].type()
                            << endl;
                    }
                }
            }
            else
            {
                Pout<< "level:" << fineLevelIndex << " : no matrix" << endl;
            }
        }
        Pout<< endl;
    }


    if ( ! matrixLevels_.size())
    {
        FatalErrorIn
        (
            "GAMGSolver::GAMGSolver"
            "("
            "const word& fieldName,"
            "const lduMatrix& matrix,"
            "const FieldField<Field, scalar>& interfaceBouCoeffs,"
            "const FieldField<Field, scalar>& interfaceIntCoeffs,"
            "const lduInterfaceFieldPtrsList& interfaces,"
            "const dictionary& solverControls"
            ")"
        )   << "No coarse levels created, either matrix too small for GAMG"
               " or nCellsInCoarsestLevel too large.\n"
               "    Either choose another solver of reduce "
               "nCellsInCoarsestLevel."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGSolver::~GAMGSolver()
{
    if (!cacheAgglomeration_)
    {
        delete &agglomeration_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::readControls()
{
    lduMatrix::solver::readControls();

    // we could also consider supplying defaults here too
    controlDict_.readIfPresent("cacheAgglomeration", cacheAgglomeration_);
    controlDict_.readIfPresent("nPreSweeps", nPreSweeps_);
    controlDict_.readIfPresent
    (
        "preSweepsLevelMultiplier",
        preSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPreSweeps", maxPreSweeps_);
    controlDict_.readIfPresent("nPostSweeps", nPostSweeps_);
    controlDict_.readIfPresent
    (
        "postSweepsLevelMultiplier",
        postSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPostSweeps", maxPostSweeps_);
    controlDict_.readIfPresent("nFinestSweeps", nFinestSweeps_);
    controlDict_.readIfPresent("interpolateCorrection", interpolateCorrection_);
    controlDict_.readIfPresent("scaleCorrection", scaleCorrection_);

    if (debug)
    {
        Pout<< "GAMGSolver settings :"
            << " cacheAgglomeration:" << cacheAgglomeration_
            << " nPreSweeps:" << nPreSweeps_
            << " preSweepsLevelMultiplier:" << preSweepsLevelMultiplier_
            << " maxPreSweeps:" << maxPreSweeps_
            << " nPostSweeps:" << nPostSweeps_
            << " postSweepsLevelMultiplier:" << postSweepsLevelMultiplier_
            << " maxPostSweeps:" << maxPostSweeps_
            << " nFinestSweeps:" << nFinestSweeps_
            << " interpolateCorrection:" << interpolateCorrection_
            << " scaleCorrection:" << scaleCorrection_
            << endl;
    }
}


const Foam::lduMatrix& Foam::GAMGSolver::matrixLevel(const label i) const
{
    if (i == 0)
    {
        return matrix_;
    }
    else
    {
        return matrixLevels_[i - 1];
    }
}


const Foam::lduInterfaceFieldPtrsList& Foam::GAMGSolver::interfaceLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaces_;
    }
    else
    {
        return interfaceLevels_[i - 1];
    }
}


const Foam::FieldField<Foam::gpuField, Foam::scalar>&
Foam::GAMGSolver::interfaceBouCoeffsLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaceBouCoeffs_;
    }
    else
    {
        return interfaceLevelsBouCoeffs_[i - 1];
    }
}


const Foam::FieldField<Foam::gpuField, Foam::scalar>&
Foam::GAMGSolver::interfaceIntCoeffsLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaceIntCoeffs_;
    }
    else
    {
        return interfaceLevelsIntCoeffs_[i - 1];
    }
}


// ************************************************************************* //
