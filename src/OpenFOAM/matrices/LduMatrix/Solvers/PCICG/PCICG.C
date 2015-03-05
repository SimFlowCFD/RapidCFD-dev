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

#include "PCICG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::PCICG<Type, DType, LUType>::PCICG
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverDict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
typename Foam::SolverPerformance<Type>
Foam::PCICG<Type, DType, LUType>::solve(gpuField<Type>& psi) const
{
    word preconditionerName(this->controlDict_.lookup("preconditioner"));

    if(preconditionerName != "none")
        preconditionerName = "diagonal";

    // --- Setup class containing solver performance data
    SolverPerformance<Type> solverPerf
    (
        preconditionerName + typeName,
        this->fieldName_
    );

    register label nCells = psi.size();


    gpuField<Type> pA(nCells);

    gpuField<Type> wA(nCells);

    Type wArA = solverPerf.great_*pTraits<Type>::one;
    Type wArAold = wArA;

    // --- Calculate A.psi
    this->matrix_.Amul(wA, psi);

    // --- Calculate initial residual field
    gpuField<Type> rA(this->matrix_.source() - wA);

    // --- Calculate normalisation factor
    Type normFactor = this->normFactor(psi, wA, pA);

    if (LduMatrix<Type, DType, LUType>::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = cmptDivide(gSumCmptMag(rA), normFactor);
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        this->minIter_ > 0
     || !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        preconPtr = LduMatrix<Type, DType, LUType>::preconditioner::New
        (
            *this,
            this->controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            preconPtr->precondition(wA, rA);

            // --- Update search directions:
            wArA = gSumCmptProd(wA, rA, matrix().mesh().comm());

            if (solverPerf.nIterations() == 0)
            {
                thrust::copy(wA.begin(),wA.end(),pA.begin());
            }
            else
            {
                Type beta = cmptDivide
                (
                    wArA,
                    stabilise(wArAold, solverPerf.vsmall_)
                );

                thrust::transform
                (
                   wA.begin(),
                   wA.end(),
                   thrust::make_transform_iterator
                   (
                       pA.begin(),
                       cmptMultiplyBinaryFunctionSFFunctor<Type,Type,Type>(beta)
                   ),
                   pA.begin(),
                   addOperatorFunctor<Type,Type,Type>()
                );
            }


            // --- Update preconditioned residual
            this->matrix_.Amul(wA, pA);

            Type wApA = gSumCmptProd(wA, pA, matrix().mesh().comm());


            // --- Test for singularity
            if
            (
                solverPerf.checkSingularity
                (
                    cmptDivide(cmptMag(wApA), normFactor)
                )
            )
            {
                break;
            }


            // --- Update solution and residual:

            Type alpha = cmptDivide
            (
                wArA,
                stabilise(wApA, solverPerf.vsmall_)
            );

            thrust::transform
            (
                psi.begin(),
                psi.end(),
                thrust::make_transform_iterator
                (
                    pA.begin(),
                    cmptMultiplyBinaryFunctionSFFunctor<Type,Type,Type>(alpha)
                ),
                psi.begin(),
                addOperatorFunctor<Type,Type,Type>()
            );

            thrust::transform
            (
                rA.begin(),
                rA.end(),
                thrust::make_transform_iterator
                (
                    wA.begin(),
                    cmptMultiplyBinaryFunctionSFFunctor<Type,Type,Type>(alpha)
                ),
                rA.begin(),
                subtractOperatorFunctor<Type,Type,Type>()
            );

            solverPerf.finalResidual() =
                cmptDivide(gSumCmptMag(rA, matrix().mesh().comm()), normFactor);

        } while
        (
            (
                solverPerf.nIterations()++ < this->maxIter_
            && !solverPerf.checkConvergence(this->tolerance_, this->relTol_)
            )
         || solverPerf.nIterations() < this->minIter_
        );
    }

    return solverPerf;
}


// ************************************************************************* //
