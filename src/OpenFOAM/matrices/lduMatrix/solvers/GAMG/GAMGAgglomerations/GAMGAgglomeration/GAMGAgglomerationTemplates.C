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

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"
#include "GAMGAgglomerationF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const labelgpuList& fineToCoarse
) const
{
    cf = pTraits<Type>::zero;

    thrust::for_each
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            fineToCoarse.begin(),
            ff.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            fineToCoarse.end(),
            ff.end()
        )),
        GAMGAgglomerationRestrictFunctor<Type>
        (
            cf.data()
        )
    );
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelgpuList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    restrictField
    (
        cf, 
        ff,
        fineToCoarse
    );
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelgpuList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictFaceField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = pTraits<Type>::zero;

    thrust::for_each
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            fineToCoarse.begin(),
            ff.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            fineToCoarse.end(),
            ff.end()
        )),
        GAMGAgglomerationRestrictFaceFunctor<Type>
        (
            cf.data()
        )
    );
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressingHost_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictFaceField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = pTraits<Type>::zero;

    forAll(fineToCoarse, ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
            cf[cFace] += ff[ffacei];
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    gpuField<Type>& ff,
    const gpuField<Type>& cf,
    const label levelIndex
) const
{
    const labelgpuList& fineToCoarse = restrictAddressing_[levelIndex];

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            cf.begin(),
            fineToCoarse.begin()
        ),
        thrust::make_permutation_iterator
        (
            cf.begin(),
            fineToCoarse.end()
        ),
        ff.begin()
    );
}


// ************************************************************************* //
