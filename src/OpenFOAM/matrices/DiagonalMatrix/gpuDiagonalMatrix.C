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

#include "gpuDiagonalMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
struct DiagonalMatrixInvertFunctor : public std::unary_function<Type,Type>
{  
    __HOST____DEVICE__
    Type operator ()(const Type& x)
    {
        if (mag(x) < VSMALL)
        {
            return Type(0);
        }
        else
        {
            return Type(1)/x;
        }
    }
};	
   
}

template<class Type>
inline Foam::gpuDiagonalMatrix<Type>::gpuDiagonalMatrix()
:
    gpuList<Type>()
{}


template<class Type>
template<class Form>
Foam::gpuDiagonalMatrix<Type>::gpuDiagonalMatrix(const Matrix<Form, Type>& a)
:
    gpuList<Type>(min(a.n(), a.m()))
{
	List<Type> tmp(this->size());
	
    forAll(tmp, i)
    {
        tmp.operator[](i) = a[i][i];
    }
    this->operator()(tmp);
}


template<class Type>
Foam::gpuDiagonalMatrix<Type>::gpuDiagonalMatrix(const label size)
:
    gpuList<Type>(size)
{}


template<class Type>
Foam::gpuDiagonalMatrix<Type>::gpuDiagonalMatrix(const label size, const Type& val)
:
    gpuList<Type>(size, val)
{}


template<class Type>
Foam::gpuDiagonalMatrix<Type>& Foam::gpuDiagonalMatrix<Type>::invert()
{
    inv(*this);

    return this;
}


template<class Type>
Foam::gpuDiagonalMatrix<Type> Foam::inv(const gpuDiagonalMatrix<Type>& A)
{
    gpuDiagonalMatrix<Type> Ainv = A;

    thrust::transform
    (
        A.begin(),
        A.end(),
        Ainv(),
        DiagonalMatrixInvertFunctor<Type>()
    );

    return Ainv;
}


// ************************************************************************* //
