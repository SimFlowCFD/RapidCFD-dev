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

#include "adjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{
template<bool positive>
struct adjustPhiFilter : public std::unary_function<scalar,scalar>{
    __HOST____DEVICE__
    scalar operator()(const scalar& s){
        if(positive)
            return s>0?s:0;
        else
            return s<0?s:0;
    } 
};

struct adjustPhiCorrectFunctor : public std::unary_function<scalar,scalar>{
    const scalar massCorr;
    adjustPhiCorrectFunctor(const scalar _massCorr):massCorr(_massCorr){}
    __HOST____DEVICE__
    scalar operator()(const scalar& phi){
        if (phi > 0.0)
            return phi * massCorr;
        else
            return phi;
    }
};
}

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p
)
{
    if (p.needReference())
    {
        // p coefficients should not be updated here
        // p.boundaryField().updateCoeffs();

        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        surfaceScalarField::GeometricBoundaryField& bphi = phi.boundaryField();

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!isA<processorFvsPatchScalarField>(phip))
            {
                massIn = thrust::reduce(thrust::make_transform_iterator(phip.begin(),adjustPhiFilter<false>()),
                                            thrust::make_transform_iterator(phip.end(),adjustPhiFilter<false>()),
                                            0,
                                            thrust::minus<scalar>());
                if (Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up))
                {
                    
                    fixedMassOut = thrust::reduce(thrust::make_transform_iterator(phip.begin(),adjustPhiFilter<true>()),
                                            thrust::make_transform_iterator(phip.end(),adjustPhiFilter<true>()),
                                            0,
                                            thrust::plus<scalar>());
                    /*
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                    */
                }
                else
                {
                    
                    adjustableMassOut = thrust::reduce(thrust::make_transform_iterator(phip.begin(),adjustPhiFilter<true>()),
                                            thrust::make_transform_iterator(phip.end(),adjustPhiFilter<true>()),
                                            0,
                                            thrust::plus<scalar>());
                    /*
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                    */
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = VSMALL + sum(mag(phi)).value();

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > VSMALL
         && magAdjustableMassOut/totalFlux > SMALL
        )
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            FatalErrorIn
            (
                "adjustPhi"
                "("
                    "surfaceScalarField&, "
                    "const volVectorField&,"
                    "volScalarField&"
                ")"
            )   << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!isA<processorFvsPatchScalarField>(phip))
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    thrust::transform(phip.begin(),phip.end(),phip.begin(),
                                      adjustPhiCorrectFunctor(massCorr));
/*
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
*/
                }
            }
        }

        return mag(massIn)/totalFlux < SMALL
            && mag(fixedMassOut)/totalFlux < SMALL
            && mag(adjustableMassOut)/totalFlux < SMALL;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
