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
#include "vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

struct multiplyTupleFunctor : public std::unary_function<thrust::tuple<scalar,scalar>,scalar>
{
    __HOST____DEVICE__
    scalar operator()(const thrust::tuple<scalar,scalar>& t)
    {
        return thrust::get<0>(t)*thrust::get<1>(t);
    }
};

struct GAMGSolverScaleFunctor
{
    const scalar sf;
    GAMGSolverScaleFunctor(const scalar sf_):sf(sf_){}
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& field, const Tuple& t)
    {
        return sf*field + (thrust::get<0>(t) - sf*thrust::get<1>(t))/thrust::get<2>(t);
    }
};

}

void Foam::GAMGSolver::scale
(
    scalargpuField& field,
    scalargpuField& Acf,
    const lduMatrix& A,
    const FieldField<gpuField, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalargpuField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );

    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    scalingFactorNum  = 
        thrust::reduce
        (
            thrust::make_transform_iterator
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    source.begin(),
                    field.begin()
                )),
                multiplyTupleFunctor()
            ),
            thrust::make_transform_iterator
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    source.end(),
                    field.end()
                )),
                multiplyTupleFunctor()
            ),
            0.0,
            thrust::plus<scalar>()
        );

    scalingFactorDenom  = 
        thrust::reduce
        (
            thrust::make_transform_iterator
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    Acf.begin(),
                    field.begin()
                )),
                multiplyTupleFunctor()
            ),
            thrust::make_transform_iterator
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    Acf.end(),
                    field.end()
                )),
                multiplyTupleFunctor()
            ),
            0.0,
            thrust::plus<scalar>()
        );

/*
    forAll(field, i)
    {
        scalingFactorNum += source[i]*field[i];
        scalingFactorDenom += Acf[i]*field[i];
    }
*/
    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    A.mesh().reduce(scalingVector, sumOp<vector2D>());

    scalar sf = scalingVector.x()/stabilise(scalingVector.y(), VSMALL);

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalargpuField& D = A.diag();

/*
    forAll(field, i)
    {
        field[i] = sf*field[i] + (source[i] - sf*Acf[i])/D[i];
    }
*/

    thrust::transform
    (
        field.begin(),
        field.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            source.begin(),
            Acf.begin(),
            D.begin()
        )),
        field.begin(),
        GAMGSolverScaleFunctor(sf)
    );
}


// ************************************************************************* //
