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

#include "jumpCyclicFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchFieldsTypeName(jumpCyclic);

    struct jumpCyclicCoeffsFunctor{
        __HOST____DEVICE__
        scalar operator()(const scalar& r, const thrust::tuple<scalar,scalar>& t){
            return r - thrust::get<0>(t)*thrust::get<1>(t);
        }
    };
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::jumpCyclicFvPatchField<Foam::scalar>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField& psiInternal,
    const scalargpuField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    scalargpuField pnf(this->size());

    const labelgpuList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    // only apply jump to original field
    if (&psiInternal == &this->internalField())
    {
        gpuField<scalar> jf(this->jump());

        if (!this->cyclicPatch().owner())
        {
            jf *= -1.0;
        }
/*
        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]] - jf[facei];
        }
*/
        thrust::transform(thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()),
                          thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.end()),
                          jf.begin(),
                          pnf.begin(),
                          subtractOperatorFunctor<scalar,scalar,scalar>());
    }
    else
    {
        thrust::copy(thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()),
                     thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.end()),
                     pnf.begin());
/*        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]];
        }*/
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    const labelgpuList& faceCells = this->cyclicPatch().faceCells();
/*    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }*/

    thrust::transform(thrust::make_permutation_iterator(result.begin(),faceCells.begin()),
                      thrust::make_permutation_iterator(result.begin(),faceCells.end()),
                      thrust::make_zip_iterator(thrust::make_tuple(coeffs.begin(),pnf.begin())),
                      thrust::make_permutation_iterator(result.begin(),faceCells.begin()),
                      jumpCyclicCoeffsFunctor());
}


// ************************************************************************* //
