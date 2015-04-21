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

#include "PCG.H"
#include "lduMatrixSolverFunctors.H"
#include "PCGCache.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCG>
        addPCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCG::PCG
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::solverPerformance Foam::PCG::solve
(
    scalargpuField& psi,
    const scalargpuField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    register label nCells = psi.size();


    scalargpuField pA(PCGCache::pA(matrix_.level(),nCells),nCells);
    scalargpuField wA(PCGCache::wA(matrix_.level(),nCells),nCells);

    scalar wArA = solverPerf.great_;
    scalar wArAold = wArA;

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalargpuField rA(PCGCache::rA(matrix_.level(),nCells),nCells);
    thrust::transform
    (
        source.begin(),
        source.end(),
        wA.begin(),
        rA.begin(),
        minusOp<scalar>()
    );

    // --- Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source, wA, pA);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = gSumMag(rA, matrix().mesh().comm())/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr =
        lduMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            preconPtr->precondition(wA, rA, cmpt);
            
            // --- Update search directions:
            wArA = gSumProd(wA, rA, matrix().mesh().comm());

            if (solverPerf.nIterations() == 0)
            {
                thrust::copy(wA.begin(),wA.end(),pA.begin());
            }
            else
            {
                scalar beta = wArA/wArAold;

                thrust::transform
                (
                    wA.begin(),
                    wA.end(),
                    pA.begin(),
                    pA.begin(),
                    wAPlusBetaPAFunctor(beta)
                );
            }


            // --- Update preconditioned residual
            matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

            scalar wApA = gSumProd(wA, pA, matrix().mesh().comm());


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and residual:

            scalar alpha = wArA/wApA;

            thrust::transform
            (
                psi.begin(),
                psi.end(),
                pA.begin(),
                psi.begin(),
                psiPlusAlphaPAFunctor(alpha)
            );

            thrust::transform
            (
                rA.begin(),
                rA.end(),
                wA.begin(),
                rA.begin(),
                rAMinusAlphaWAFunctor(alpha)
            );

            solverPerf.finalResidual() = gSumMag(rA, matrix().mesh().comm())/normFactor;

        } while
        (
            (
                solverPerf.nIterations()++ < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    return solverPerf;
}


// ************************************************************************* //
