/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "PBiCICG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::PBiCICG<Type, DType, LUType>::PBiCICG
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
Foam::SolverPerformance<Type>
Foam::PBiCICG<Type, DType, LUType>::solve(gpuField<Type>& psi) const
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

    gpuField<Type> pT(nCells, pTraits<Type>::zero);

    gpuField<Type> wA(nCells);

    gpuField<Type> wT(nCells);

    Type wArT = solverPerf.great_*pTraits<Type>::one;
    Type wArTold = wArT;

    // --- Calculate A.psi and T.psi
    this->matrix_.Amul(wA, psi);
    this->matrix_.Tmul(wT, psi);

    // --- Calculate initial residual and transpose residual fields
    gpuField<Type> rA(this->matrix_.source() - wA);
    gpuField<Type> rT(this->matrix_.source() - wT);

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
    if (!solverPerf.checkConvergence(this->tolerance_, this->relTol_))
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
            // --- Store previous wArT
            wArTold = wArT;

            // --- Precondition residuals
            preconPtr->precondition(wA, rA);
            preconPtr->preconditionT(wT, rT);

            // --- Update search directions:
            wArT = gSumCmptProd(wA, rT);

            if (solverPerf.nIterations() == 0)
            {
                thrust::copy(wA.begin(),wA.end(),pA.begin());
                thrust::copy(wT.begin(),wT.end(),pT.begin());
            }
            else
            {
                Type beta = cmptDivide
                (
                    wArT,
                    stabilise(wArTold, solverPerf.vsmall_)
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

               thrust::transform
               (
                   wT.begin(),
                   wT.end(),
                   thrust::make_transform_iterator
                   (
                       pT.begin(),
                       cmptMultiplyBinaryFunctionSFFunctor<Type,Type,Type>(beta)
                   ),
                   pT.begin(),
                   addOperatorFunctor<Type,Type,Type>()
               );
            }


            // --- Update preconditioned residuals
            this->matrix_.Amul(wA, pA);
            this->matrix_.Tmul(wT, pT);

            Type wApT = gSumCmptProd(wA, pT);

            // --- Test for singularity
            if
            (
                solverPerf.checkSingularity
                (
                    cmptDivide(cmptMag(wApT), normFactor)
                )
            )
            {
                break;
            }


            // --- Update solution and residual:

            Type alpha = cmptDivide
            (
                wArT,
                stabilise(wApT, solverPerf.vsmall_)
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

            thrust::transform
            (
                rT.begin(),
                rT.end(),
                thrust::make_transform_iterator
                (
                    wT.begin(),
                    cmptMultiplyBinaryFunctionSFFunctor<Type,Type,Type>(alpha)
                ),
                rT.begin(),
                subtractOperatorFunctor<Type,Type,Type>()
            );

            solverPerf.finalResidual() =
                cmptDivide(gSumCmptMag(rA), normFactor);

        } while
        (
            solverPerf.nIterations()++ < this->maxIter_
        && !(solverPerf.checkConvergence(this->tolerance_, this->relTol_))
        );
    }

    return solverPerf;
}


// ************************************************************************* //
