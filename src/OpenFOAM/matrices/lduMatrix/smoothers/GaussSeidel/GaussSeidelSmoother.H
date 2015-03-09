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

Class
    Foam::GaussSeidelSmoother

Description
    Dummy Gauss-Seidel as a proxy to Jacobi

SourceFiles
    GaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef GaussSeidelSmoother_H
#define GaussSeidelSmoother_H

#include "lduMatrix.H"
#include "JacobiSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class GaussSeidelSmoother
:
    public JacobiSmoother
{

public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from components
        GaussSeidelSmoother
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<gpuField, scalar>& interfaceBouCoeffs,
            const FieldField<gpuField, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
