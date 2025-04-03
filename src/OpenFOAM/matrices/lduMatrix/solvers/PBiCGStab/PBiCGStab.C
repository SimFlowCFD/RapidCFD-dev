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

#include "PBiCGStab.H"
#include "lduMatrixSolverFunctors.H"
#include "PCGCache.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PBiCGStab, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<PBiCGStab>
        addPBiCGStabAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBiCGStab::PBiCGStab
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

Foam::solverPerformance Foam::PBiCGStab::solve
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

    scalargpuField yA(PCGCache::yA(matrix_.level(),nCells),nCells);

    // --- Calculate A.psi
    matrix_.Amul(yA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalargpuField rA(PCGCache::rA(matrix_.level(),nCells),nCells);

    thrust::transform
    (
        source.begin(),
        source.end(),
        yA.begin(),
        rA.begin(),
        minusOp<scalar>()
    );

    // --- Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source, yA, pA);

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
        scalargpuField AyA(PCGCache::AyA(matrix_.level(),nCells),nCells);
        scalargpuField sA(PCGCache::sA(matrix_.level(),nCells),nCells);
        scalargpuField zA(PCGCache::zA(matrix_.level(),nCells),nCells);
        scalargpuField tA(PCGCache::tA(matrix_.level(),nCells),nCells);
        scalargpuField result1(PCGCache::result1(matrix_.level(),nCells),nCells);
 
        // --- Store initial residual
        const scalargpuField rA0(rA);
    
        // --- Initial values not used
        scalar rA0rA = 0;
        scalar alpha = 0;
        scalar omega = 0;

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
            // --- Store previous rA0rA
            const scalar rA0rAold = rA0rA;

            rA0rA = gSumProd(rA0, rA, matrix().mesh().comm());

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(rA0rA)))
            {
                break;
            }
            // --- Update pA
            if (solverPerf.nIterations() == 0)
            {
                thrust::copy(rA.begin(),rA.end(),pA.begin());
            }
            else
            {
                // --- Test for singularity
                if (solverPerf.checkSingularity(mag(omega)))
                {
                    break;
                }

                const scalar beta = (rA0rA/rA0rAold)*(alpha/omega);

// pAPtr[cell] = rAPtr[cell] + beta*(pAPtr[cell] - omega*AyAPtr[cell]);
// Split this into:
//  result1 = pAPtr[cell] - omega*AyAPtr[cell]
//  pAPtr[cell] = rAPtr[cell] + beta*result1
                thrust::transform
                (
                    pA.begin(),
                    pA.end(),
                    AyA.begin(),
                    result1.begin(),
                    pAMinusOmegaAYAFunctor(omega)
                );

                thrust::transform
                (
                    rA.begin(),
                    rA.end(),
                    result1.begin(),
                    pA.begin(),
                    rAPlusBetaResult1Functor(beta)
                );
            }

            // --- Precondition pA
            preconPtr->precondition(yA, pA, cmpt);

            // --- Calculate AyA
            matrix_.Amul(AyA, yA, interfaceBouCoeffs_, interfaces_, cmpt);

            const scalar rA0AyA = gSumProd(rA0, AyA, matrix().mesh().comm());

            alpha = rA0rA/rA0AyA;

            // --- Calculate sA
                thrust::transform
                (
                    rA.begin(),
                    rA.end(),
                    AyA.begin(),
                    sA.begin(),
                    rAMinusAlphaAYAFunctor(alpha)
                );

            // --- Test sA for convergence
            solverPerf.finalResidual() =
                gSumMag(sA, matrix().mesh().comm())/normFactor;


            if (solverPerf.checkConvergence(tolerance_, relTol_))
            {
                thrust::transform
                (
                    psi.begin(),
                    psi.end(),
                    yA.begin(),
                    psi.begin(),
                    psiPlusAlphaYAFunctor(alpha)
                );

                solverPerf.nIterations()++;

                return solverPerf;
            }

            // --- Precondition sA
            preconPtr->precondition(zA, sA, cmpt);

            // --- Calculate tA
            matrix_.Amul(tA, zA, interfaceBouCoeffs_, interfaces_, cmpt);

            const scalar tAtA = gSumProd(tA, tA, matrix().mesh().comm());

            // --- Calculate omega from tA and sA
            //     (cheaper than using zA with preconditioned tA)
            omega = gSumProd(tA, sA, matrix().mesh().comm())/tAtA;

            // --- Update solution and residual

//psiPtr[cell] += alpha*yAPtr[cell] + omega*zAPtr[cell];
// split into
//  psiPtr[cell] += alpha*yAPtr[cell]
//  psiPtr[cell] += omega*zAPtr[cell]

            thrust::transform
            (
                psi.begin(),
                psi.end(),
                yA.begin(),
                psi.begin(),
                psiPlusAlphaYAFunctor(alpha)
            );

            thrust::transform
            (
                psi.begin(),
                psi.end(),
                yA.begin(),
                psi.begin(),
                psiPlusOmegaZAFunctor(omega)
            );

// rAPtr[cell] = sAPtr[cell] - omega*tAPtr[cell];

            thrust::transform
            (
                sA.begin(),
                sA.end(),
                tA.begin(),
                rA.begin(),
                sAMinusOmegaTAFunctor(omega)
            );

            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;

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
