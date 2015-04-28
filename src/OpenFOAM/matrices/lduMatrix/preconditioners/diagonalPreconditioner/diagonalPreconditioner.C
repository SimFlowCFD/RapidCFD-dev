/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "diagonalPreconditioner.H"
#include "lduMatrixSolutionCache.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diagonalPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<diagonalPreconditioner>
        adddiagonalPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<diagonalPreconditioner>
        adddiagonalPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diagonalPreconditioner::diagonalPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD
    (
        lduMatrixSolutionCache::first(sol.matrix().diag().size()),
        sol.matrix().diag().size()
    )
{ 
    const scalargpuField& Diag = solver_.matrix().diag();

    thrust::transform
    (
        Diag.begin(),
        Diag.end(),
        rD.begin(),
        divideOperatorSFFunctor<scalar,scalar,scalar>(1.0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diagonalPreconditioner::precondition
(
    scalargpuField& wA,
    const scalargpuField& rA,
    const direction
) const
{
    thrust::transform
    (
        rD.begin(),
        rD.end(),
        rA.begin(),
        wA.begin(),
        multiplyOperatorFunctor<scalar,scalar,scalar>()
    );
}


// ************************************************************************* //
