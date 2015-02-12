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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::gatherList
(
    const label comm,
    const labelList& procIDs,

    const Type& myVal,
    List<Type>& allVals,
    const int tag
)
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allVals.setSize(procIDs.size());

        allVals[0] = myVal;
        for (label i = 1; i < procIDs.size(); i++)
        {
            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,
                tag,
                comm
            );

            fromSlave >> allVals[i];
        }
    }
    else
    {
        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            tag,
            comm
        );
        toMaster << myVal;
    }
}

namespace Foam
{

template<class Type>
struct GAMGAgglomerationRestrictFunctor
{
    const Type* ff;
    const label* sort;
    const Type zero;

    GAMGAgglomerationRestrictFunctor
    (
        const Type* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort),
        zero(pTraits<Type>::zero)
    {}

    __HOST____DEVICE__
    Type operator()(const label& start, const label& end)
    {
        Type out = zero;

        for(label i = start; i<end; i++)
        {
            out += ff[sort[i]];
        }

        return out;
    }

};

}

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
        GAMGAgglomerationRestrictFunctor<Type>
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
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelgpuList& sort = restrictSortAddressing_[fineLevelIndex];
    const labelgpuList& target = restrictTargetAddressing_[fineLevelIndex];
    const labelgpuList& targetStart = restrictTargetStartAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != sort.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictField"
            "(Field<Type>& cf, const Field<Type>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << sort.size()
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

    label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label fineComm = UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<int>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        Field<Type> cfh(cf.size());

        globalIndex::gather
        (
            offsets,
            fineComm,
            procIDs,
            cfh,
            UPstream::msgType(),
            Pstream::nonBlocking    //Pstream::scheduled
        );

        thrust::copy
        (
            cfh.begin(),
            cfh.end(),
            cf.begin()
        );
    }
}

namespace Foam
{

template<class Type>
struct nonNegativeGAMGFunctor
{
    __HOST____DEVICE__
    bool operator()(const Type& x)
    {
        return x >= 0;
    }
};

}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    gpuField<Type>& cf,
    const gpuField<Type>& ff,
    const label fineLevelIndex
) const
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
            << " level = " << sort.size()
            << abort(FatalError);
    }

    cf = pTraits<Type>::zero;

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
        GAMGAgglomerationRestrictFunctor<Type>
        (
            ff.data(),
            sort.data()
        ),
        nonNegativeGAMGFunctor<label>()
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
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelgpuList& fineToCoarse = restrictAddressing_[levelIndex];

    label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label coarseComm = UPstream::parent
        (
            procCommunicator_[coarseLevelIndex]
        );

        const List<int>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        label localSize = nCells_[levelIndex];

        Field<Type> allCf(localSize);
        Field<Type> cfh(cf.size());
        
        thrust::copy(cf.begin(),cf.end(),cfh.begin());

        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cfh,
            allCf,
            UPstream::msgType(),
            Pstream::nonBlocking    //Pstream::scheduled
        );

        gpuField<Type> allCfGpu(allCf);

        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                allCfGpu.begin(),
                fineToCoarse.begin()
            ),
            thrust::make_permutation_iterator
            (
                allCfGpu.begin(),
                fineToCoarse.end()
            ),
            ff.begin()
        );
    }
    else
    {
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
}


// ************************************************************************* //
