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
    Foam::DiagonalMatrix

Description
    gpuDiagonalMatrix<Type> is a 2D diagonal matrix of objects
    of type Type, size nxn

SourceFiles
    gpuDiagonalMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef gpuDiagonalMatrix_H
#define gpuDiagonalMatrix_H

#include "gpuList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * *  * * * * * * Class Forward declaration  * * * * * * * * * * * //

template<class Form, class Type> class Matrix;

/*---------------------------------------------------------------------------*\
                           Class DiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gpuDiagonalMatrix
:
    public gpuList<Type>
{
public:

    // Constructors

        //- Null constructor.
        gpuDiagonalMatrix<Type>();

        //- Construct from diagonal component of a Matrix
        template<class Form>
        gpuDiagonalMatrix<Type>(const Matrix<Form, Type>&);

        //- Construct empty from size
        gpuDiagonalMatrix<Type>(const label size);

        //- Construct from size and a value
        gpuDiagonalMatrix<Type>(const label, const Type&);


    // Member functions

        //- Invert the diagonal matrix and return itself
        gpuDiagonalMatrix<Type>& invert();
};


// Global functions

//- Return the diagonal Matrix inverse
template<class Type>
gpuDiagonalMatrix<Type> inv(const gpuDiagonalMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gpuDiagonalMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
