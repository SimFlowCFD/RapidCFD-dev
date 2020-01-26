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
#include "DeviceConfig.H"
#include "GAMGAgglomerationF.H"
#include "GAMGAgglomerateF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const labelgpuList& sort,
    const labelgpuList& target,
    const labelgpuList& targetStart
) const
{
    cf = pTraits<Type>::zero;

    thrust::transform
    (
        targetStart.begin(),
        targetStart.end()-1,
        targetStart.begin()+1,
        thrust::make_permutation_iterator
        (
            cf.begin(),
            target.begin()
        ),
        GAMG::restrict<Type>
        (
            ff.data(),
            sort.data()
        )
    );
}

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const labelgpuList& fineToCoarse
) const
{
    if(!is_number<Type>::value || !useAtomic())
    {
        FatalErrorIn
            (
                "void GAMGAgglomeration::restrictField"
                "(gpuField<Type>& cf, const gpuField<Type>& ff, "
                "const labelgpuList& fineToCoarse) const"
            )   << "atempting to perform atomic restriction"
                << " when it is not supported " 
                << abort(FatalError);
    }
   
    // TODO maybe it is redundant?
    cf = pTraits<Type>::zero;
    
    auto restrictOp = GAMG::atomicRestrict<Type>(
        cf.data(),
        ff.data(),
        fineToCoarse.data()
    );
    auto start = thrust::make_counting_iterator(0);

    thrust::for_each(start, start + fineToCoarse.size(), restrictOp);
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const label fineLevelIndex
) const
{
    if (is_number<Type>::value && useAtomic())
    {
        const labelgpuList& fineToCoarse = restrictAddressing_[fineLevelIndex];
        if (ff.size() != fineToCoarse.size())
        {
            FatalErrorIn
            (
                "void GAMGAgglomeration::restrictField"
                "(gpuField<Type>& cf, const gpuField<Type>& ff, "
                "const label fineLevelIndex) const"
            )   << "field does not correspond to level " << fineLevelIndex
                << " sizes: field = " << ff.size()
                << " level = " << fineToCoarse.size()
                << abort(FatalError);
        }

        restrictField(cf, ff, fineToCoarse);
    }
    else
    {
        const labelgpuList& sort = restrictSortAddressing_[fineLevelIndex];
        const labelgpuList& target = restrictTargetAddressing_[fineLevelIndex];
        const labelgpuList& targetStart = restrictTargetStartAddressing_[fineLevelIndex];

        if (ff.size() != sort.size())
        {
            FatalErrorIn
            (
                "void GAMGAgglomeration::restrictField"
                "(Field<Type>& cf, const Field<Type>& ff, "
                "const label fineLevelIndex) const"
            )   << "field does not correspond to level " << fineLevelIndex
                << " sizes: field = " << ff.size()
                << " sort = " << sort.size()
                << abort(FatalError);
        }

        restrictField
        (
            cf, 
            ff,
            sort,
            target,
            targetStart
        );
    }
}

template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const label fineLevelIndex
) const
{
    cf = pTraits<Type>::zero;

    if (is_number<Type>::value && useAtomic())
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

        Type* cfPtr = cf.data();
        const Type* ffPtr = ff.data();

        auto restrictFace = GAMG::atomicFaceRestrict<Type>(
            cf.data(),
            ff.data(),
            fineToCoarse.data()
        );

        auto fineFaceStart = thrust::make_counting_iterator(0);
        auto fineFaceEnd = fineFaceStart + fineToCoarse.size();

        thrust::for_each(fineFaceStart, fineFaceEnd, restrictFace);
    }
    else
    {
        const labelgpuList& sort = faceRestrictSortAddressing_[fineLevelIndex];
        const labelgpuList& target = faceRestrictTargetAddressing_[fineLevelIndex];
        const labelgpuList& targetStart = faceRestrictTargetStartAddressing_[fineLevelIndex];

        if (ff.size() != sort.size())
        {
            FatalErrorIn
            (
                "void GAMGAgglomeration::restrictFaceField"
                "(Field<Type>& cf, const Field<Type>& ff, "
                "const label fineLevelIndex) const"
            )   << "field does not correspond to level " << fineLevelIndex
                << " sizes: field = " << ff.size()
                << " sort = " << sort.size()
                << abort(FatalError);
        }

        thrust::transform
        (
            targetStart.begin(),
            targetStart.end()-1,
            targetStart.begin()+1,
            target.begin(),
            thrust::make_permutation_iterator
            (
                cf.begin(),
                target.begin()
            ),
            GAMG::restrict<Type>
            (
                ff.data(),
                sort.data()
            ),
            GAMG::nonNegative()
        );
    }
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
    if (is_number<Type>::value && useAtomic())
    {
        const labelgpuList& fineToCoarse = restrictAddressing_[levelIndex];
        auto coarse = thrust::make_permutation_iterator(cf.begin(), fineToCoarse.begin());

        thrust::copy(coarse, coarse+fineToCoarse.size(), ff.begin());
    }
    else
    {
        const labelgpuList& sort = restrictSortAddressing_[levelIndex];
        const labelgpuList& target = restrictTargetAddressing_[levelIndex];
        const labelgpuList& targetStart = restrictTargetStartAddressing_[levelIndex];

        thrust::for_each
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+target.size(),
            GAMG::prolong<Type>
            (
                ff.data(),
                cf.data(),
                sort.data(),
                target.data(),
                targetStart.data()
            )
        );
    }
}


// ************************************************************************* //
