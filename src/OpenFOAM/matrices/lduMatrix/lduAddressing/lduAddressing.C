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

#include "lduAddressing.H"
#include "demandDrivenData.H"
#include "scalarField.H"
#include "DynamicList.H"
#include "error.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduAddressing::calcPatchSort() const
{
    if (patchSortAddr_.size() == nPatches())
    {
        FatalErrorIn("lduAddressing::calcPatchSort() const")
            << "patch sort already calculated"
            << abort(FatalError);
    }

    patchSortAddr_.setSize(nPatches());
    patchSortCells_.setSize(nPatches());

    for(label i = 0; i < nPatches(); i++)
    {
        if( ! patchAvailable(i))
            continue;

        const labelgpuList& nbr = patchAddr(i);
        labelgpuList* sortPtr_ = new labelgpuList(nbr.size(), -1);
        patchSortAddr_.set(i,sortPtr_);

        labelgpuList& lst = *sortPtr_;

        labelgpuList nbrTmp(nbr);

        thrust::counting_iterator<label> first(0);
        thrust::copy
        (
            first,
            first+nbr.size(),
            lst.begin()
        );

        thrust::stable_sort_by_key
        (
            nbrTmp.begin(),
            nbrTmp.end(),
            lst.begin()
        );
         
        labelgpuList* cellsSortPtr= new labelgpuList(nbr.size());
        patchSortCells_.set(i,cellsSortPtr);
        labelgpuList& cellsSort = *cellsSortPtr;

        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                nbr.begin(),
                lst.begin()
            ),
            thrust::make_permutation_iterator
            (
                nbr.begin(),
                lst.end()
            ),
            cellsSort.begin()
        );
                      
        cellsSort.setSize
        (
            thrust::unique(cellsSort.begin(),cellsSort.end()) - cellsSort.begin()
        );
    }
}


void Foam::lduAddressing::calcPatchSortStart() const
{
    if (patchSortStartAddr_.size() == nPatches())
    {
        FatalErrorIn("lduAddressing::calcLosortStart() const")
            << "losort start already calculated"
            << abort(FatalError);
    }

    patchSortStartAddr_.setSize(nPatches());

    for(label i = 0; i < nPatches(); i++)
    {
        if( ! patchAvailable(i))
            continue;

        const labelgpuList& nbr = patchAddr(i);

        const labelgpuList& lsrt = patchSortAddr(i);

        labelgpuList* patchSortStartPtr_ = new labelgpuList(nbr.size() + 1, nbr.size());

        patchSortStartAddr_.set(i,patchSortStartPtr_);

        labelgpuList& lsrtStart = *patchSortStartPtr_;

        labelgpuList ones(nbr.size(),1);
        labelgpuList tmpCell(nbr.size());
        labelgpuList tmpSum(nbr.size());

        labelgpuList nbrSort(nbr.size());

        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                nbr.begin(),
                lsrt.begin()
            ),
            thrust::make_permutation_iterator
            (
                nbr.begin(),
                lsrt.end()
            ),
            nbrSort.begin()
        );

        thrust::reduce_by_key
        (
            nbrSort.begin(),
            nbrSort.end(),
            ones.begin(),
            tmpCell.begin(),
            tmpSum.begin()
        );

        thrust::exclusive_scan
        (
            tmpSum.begin(),
            tmpSum.end(),
            lsrtStart.begin()
        );
    }
}

void Foam::lduAddressing::calcLosort() const
{
    if (losortPtr_)
    {
        FatalErrorIn("lduAddressing::calcLosort() const")
            << "losort already calculated"
            << abort(FatalError);
    }

    const labelgpuList& nbr = upperAddr();
    losortPtr_ = new labelgpuList(nbr.size());

    labelgpuList& lst = *losortPtr_;

    labelgpuList nbrTmp(nbr);

    thrust::counting_iterator<label> first(0);
    thrust::copy
    (
        first,
        first+nbr.size(),
        lst.begin()
    );

    thrust::stable_sort_by_key
    (
        nbrTmp.begin(),
        nbrTmp.end(),
        lst.begin()
    );                
}


void Foam::lduAddressing::calcOwnerStart() const
{
    if (ownerStartPtr_)
    {
        FatalErrorIn("lduAddressing::calcOwnerStart() const")
            << "owner start already calculated"
            << abort(FatalError);
    }

    const labelgpuList& own = lowerAddr();

    ownerStartPtr_ = new labelgpuList(size() + 1, own.size());

    labelgpuList& ownStart = *ownerStartPtr_;

    labelgpuList ones(own.size()+size(),1);
    labelgpuList tmpCell(size());
    labelgpuList tmpSum(size());
    
    labelgpuList ownSort(own.size()+size());
    
    
    thrust::copy
    (
        own.begin(),
        own.end(),
        ownSort.begin()
    );
                 
    thrust::copy
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+size(),
        ownSort.begin()+own.size()
    );
                 
    thrust::fill
    (
        ones.begin()+own.size(),
        ones.end(),
        0
    );
    
    thrust::stable_sort_by_key
    (
        ownSort.begin(),
        ownSort.end(),
        ones.begin()
    );  
    

    thrust::reduce_by_key
    (
        ownSort.begin(),
        ownSort.end(),
        ones.begin(),
        tmpCell.begin(),
        tmpSum.begin()
    );

    thrust::exclusive_scan
    (
        tmpSum.begin(),
        tmpSum.end(),
        ownStart.begin()
    );
}


void Foam::lduAddressing::calcLosortStart() const
{
    if (losortStartPtr_)
    {
        FatalErrorIn("lduAddressing::calcLosortStart() const")
            << "losort start already calculated"
            << abort(FatalError);
    }

    const labelgpuList& nbr = upperAddr();

    const labelgpuList& lsrt = losortAddr();

    losortStartPtr_ = new labelgpuList(size() + 1, nbr.size());

    labelgpuList& lsrtStart = *losortStartPtr_;

    labelgpuList ones(nbr.size()+size(),1);
    labelgpuList tmpCell(size());
    labelgpuList tmpSum(size());

    labelgpuList nbrSort(nbr.size()+size());

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            nbr.begin(),
            lsrt.begin()
        ),
        thrust::make_permutation_iterator
        (
            nbr.begin(),
            lsrt.end()
        ),
        nbrSort.begin()
    );
                 
    thrust::copy
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+size(),
        nbrSort.begin()+nbr.size()
    );
                 
    thrust::fill
    (
        ones.begin()+nbr.size(),
        ones.end(),
        0
    );
                 

    thrust::stable_sort_by_key
    (
        nbrSort.begin(),
        nbrSort.end(),
        ones.begin()
    );     

    thrust::reduce_by_key
    (
        nbrSort.begin(),
        nbrSort.end(),
        ones.begin(),
        tmpCell.begin(),
        tmpSum.begin()
    );

    thrust::exclusive_scan
    (
        tmpSum.begin(),
        tmpSum.begin()+size(),
        lsrtStart.begin()
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lduAddressing::~lduAddressing()
{
    deleteDemandDrivenData(losortPtr_);
    deleteDemandDrivenData(ownerStartPtr_);
    deleteDemandDrivenData(losortStartPtr_);
    deleteDemandDrivenData(ownerSortAddrPtr_);
    
    patchSortCells_.clear();
    patchSortAddr_.clear();
    patchSortStartAddr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelgpuList& Foam::lduAddressing::losortAddr() const
{
    if (!losortPtr_)
    {
        calcLosort();
    }

    return *losortPtr_;
}

const Foam::labelgpuList& Foam::lduAddressing::ownerSortAddr() const
{
    if ( ! ownerSortAddrPtr_)
    { 
        const labelgpuList& own = lowerAddr();
        const labelgpuList& lsrt = losortAddr();

        ownerSortAddrPtr_ = new labelgpuList(own.size());
        labelgpuList& ownSort = *ownerSortAddrPtr_;
        
        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                own.begin(),
                lsrt.begin()
            ),
            thrust::make_permutation_iterator
            (
                own.begin(),
                lsrt.end()
            ),
            ownSort.begin()
        );
    }

    return *ownerSortAddrPtr_;
}

const Foam::labelgpuList& Foam::lduAddressing::ownerStartAddr() const
{
    if (!ownerStartPtr_)
    {
        calcOwnerStart();
    }

    return *ownerStartPtr_;
}


const Foam::labelgpuList& Foam::lduAddressing::losortStartAddr() const
{
    if (!losortStartPtr_)
    {
        calcLosortStart();
    }

    return *losortStartPtr_;
}

const Foam::labelgpuList& Foam::lduAddressing::patchSortCells(const label i) const
{
    if (patchSortCells_.size() != nPatches())
    {
        calcPatchSort();
    }

    return patchSortCells_[i];
}

const Foam::labelgpuList& Foam::lduAddressing::patchSortAddr(const label i) const
{
    if (patchSortAddr_.size() != nPatches())
    {
        calcPatchSort();
    }

    return patchSortAddr_[i];
}

const Foam::labelgpuList& Foam::lduAddressing::patchSortStartAddr(const label i) const
{
    if (patchSortStartAddr_.size() != nPatches())
    {
        calcPatchSortStart();
    }

    return patchSortStartAddr_[i];
}

Foam::Tuple2<Foam::label, Foam::scalar> Foam::lduAddressing::band() const
{
    const labelgpuList& owner = lowerAddr();
    const labelgpuList& neighbour = upperAddr();

    labelgpuList cellBandwidth(size(), 0);
    labelgpuList diffs(neighbour.size(),0);

    thrust::transform
    (
        neighbour.begin(),
        neighbour.end(),
        owner.begin(),
        diffs.begin(),
        subtractOperatorFunctor<label,label,label>()
    );

    thrust::transform
    (
        diffs.begin(),
        diffs.end(),
        thrust::make_permutation_iterator
        (
            cellBandwidth.begin(),
            neighbour.begin()
        ),
        thrust::make_permutation_iterator
        (
            cellBandwidth.begin(),
            neighbour.begin()
        ),
        maxBinaryFunctionFunctor<label,label,label>()
    );

    label bandwidth = max(cellBandwidth);

    // Do not use field algebra because of conversion label to scalar
    scalar profile = 
        thrust::reduce
        (
            cellBandwidth.begin(),
            cellBandwidth.end()
        );

    return Tuple2<label, scalar>(bandwidth, profile);
}


// ************************************************************************* //
