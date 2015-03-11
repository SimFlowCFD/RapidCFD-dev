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

#include "MULES.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::MULES::explicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    explicitSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}


void Foam::MULES::explicitLTSSolve
(
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const scalar psiMax,
    const scalar psiMin
)
{
    explicitLTSSolve
    (
        geometricOneField(),
        psi,
        phi,
        phiPsi,
        zeroField(), zeroField(),
        psiMax, psiMin
    );
}

namespace Foam
{

struct sumPosSumNegMULESFunctor
{
    __HOST____DEVICE__
    thrust::tuple<scalar,scalar> operator()
    (
        const scalar& phiPsiCorrs, 
        const thrust::tuple<scalar,scalar>& t
    )
    {
        scalar sumPos = thrust::get<0>(t);
        scalar sumNeg = thrust::get<1>(t);

        if (phiPsiCorrs > 0)
        {
            sumPos += phiPsiCorrs;
        }
        else
        {
            sumNeg += phiPsiCorrs;
        }

        return thrust::make_tuple(sumPos,sumNeg);
    }
};

struct lambdaMULESFunctor
{
    __HOST____DEVICE__
    thrust::tuple<scalar,bool,bool> operator()(const scalar& sumPos,const scalar& sumNeg)
    {
        scalar sum = sumPos + sumNeg;
        scalar lambda = 1;
        bool positive = false;
        bool negative = false;

        if (sum > 0 && sumPos > VSMALL)
        {
            lambda = -sumNeg/sumPos;
            positive = true;
        }
        else if (sum < 0 && sumNeg < -VSMALL)
        {
            lambda = -sumPos/sumNeg;
            negative = true;
        }

        return thrust::make_tuple(lambda,positive,negative);
    }
};

struct phiPsiCorrsMULESFunctor
{
    __HOST____DEVICE__
    scalar operator()(const scalar& phiPsiCorrs, const thrust::tuple<scalar,bool,bool>& t)
    {
        scalar out = phiPsiCorrs; 
        const scalar lambda = thrust::get<0>(t);
        const bool positive = thrust::get<1>(t);
        const bool negative = thrust::get<2>(t);

        if (positive)
        {
            if(phiPsiCorrs > 0)
            {
                out = phiPsiCorrs * lambda;
            }
        }
        else if (negative)
        {
            if(phiPsiCorrs > 0)
            {
                out = phiPsiCorrs * lambda;
            }
        }
        
        return out;
    } 
};

}

void Foam::MULES::limitSum(UPtrList<scalargpuField>& phiPsiCorrs)
{
    label size = phiPsiCorrs[0].size();

    scalargpuField sumPos(size,0);
    scalargpuField sumNeg(size,0);

    scalargpuField lambda(size,0);

    boolgpuList positive(size,false);
    boolgpuList negative(size,false);

    for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
    {
        thrust::transform
        (
            phiPsiCorrs[phasei].begin(),
            phiPsiCorrs[phasei].end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                sumPos.begin(),
                sumNeg.begin()
            )),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                sumPos.begin(),
                sumNeg.begin()
            )),
            sumPosSumNegMULESFunctor()
        );
    }


    thrust::transform
    (
        sumPos.begin(),
        sumPos.end(),
        sumNeg.begin(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            lambda.begin(),
            positive.begin(),
            negative.begin()
        )),
        lambdaMULESFunctor()
    );


    for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
    {
        thrust::transform
        (
            phiPsiCorrs[phasei].begin(),
            phiPsiCorrs[phasei].end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                lambda.begin(),
                positive.begin(),
                negative.begin()
            )),
            phiPsiCorrs[phasei].begin(),
            phiPsiCorrsMULESFunctor()
        );
    }
}


// ************************************************************************* //
