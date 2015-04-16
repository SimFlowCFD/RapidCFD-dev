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

#include "GAMGInterface.H"
#include "GAMGInterfaceF.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGInterface, 0);
    defineRunTimeSelectionTable(GAMGInterface, lduInterface);
    defineRunTimeSelectionTable(GAMGInterface, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGInterface::GAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    index_(index),
    coarseInterfaces_(coarseInterfaces),
    faceCellsHost_(is),
    faceCells_(faceCellsHost_),
    faceRestrictAddressingHost_(is),
    faceRestrictAddressing_(faceRestrictAddressingHost_)
{
    updateAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GAMGInterface::updateAddressing()
{
    GAMGAgglomeration::createSort
    (
        faceRestrictAddressing_,
        faceRestrictSortAddressing_
    );

    GAMGAgglomeration::createTarget
    (
        faceRestrictAddressing_,
        faceRestrictSortAddressing_,
        faceRestrictTargetAddressing_,
        faceRestrictTargetStartAddressing_
    );

    GAMGAgglomeration::createSort
    (
        faceCells_,
        sortCells_
    );

    GAMGAgglomeration::createTarget
    (
        faceCells_,
        sortCells_,
        cellFaces_,
        cellFacesStart_
    );
}

void Foam::GAMGInterface::combine(const GAMGInterface& coarseGi)
{
    const labelgpuList& coarseFra = coarseGi.faceRestrictAddressing_;
    const labelList& coarseFraHost = coarseGi.faceRestrictAddressingHost_;
  
    labelgpuList restTmp(faceRestrictAddressing_);

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            coarseFra.begin(),
            restTmp.begin()
        ),
        thrust::make_permutation_iterator
        (
            coarseFra.begin(),
            restTmp.end()
        ),
        faceRestrictAddressing_.begin()
    );

    forAll(faceRestrictAddressingHost_, ffi)
    {
        faceRestrictAddressingHost_[ffi] = coarseFraHost[faceRestrictAddressingHost_[ffi]];
    }

    faceCells_ = coarseGi.faceCells_;
    faceCellsHost_ = coarseGi.faceCellsHost_;

    updateAddressing();
}


Foam::tmp<Foam::labelField> Foam::GAMGInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    tmp<labelField> tresult(new labelField(size()));
    labelField& result = tresult();

    forAll(result, elemI)
    {
        result[elemI] = internalData[faceCellsHost_[elemI]];
    }

    return tresult;
}


Foam::tmp<Foam::scalargpuField> Foam::GAMGInterface::agglomerateCoeffs
(
    const scalargpuField& fineCoeffs
) const
{
    tmp<scalargpuField> tcoarseCoeffs(new scalargpuField(size(), 0.0));
    scalargpuField& coarseCoeffs = tcoarseCoeffs();

    if (fineCoeffs.size() != faceRestrictAddressing_.size())
    {
        FatalErrorIn
        (
            "GAMGInterface::agglomerateCoeffs(const scalarField&) const"
        )   << "Size of coefficients " << fineCoeffs.size()
            << " does not correspond to the size of the restriction "
            << faceRestrictAddressing_.size()
            << abort(FatalError);
    }
    if (debug && max(faceRestrictAddressing_) > size())
    {
        FatalErrorIn
        (
            "GAMGInterface::agglomerateCoeffs(const scalargpuField&) const"
        )   << "Face restrict addressing addresses outside of coarse interface"
            << " size. Max addressing:" << max(faceRestrictAddressing_)
            << " coarse size:" << size()
            << abort(FatalError);
    }

    thrust::transform
    (
        faceRestrictTargetStartAddressing_.begin(),
        faceRestrictTargetStartAddressing_.end()-1,
        faceRestrictTargetStartAddressing_.begin()+1,
        thrust::make_permutation_iterator
        (
            coarseCoeffs.begin(),
            faceRestrictTargetAddressing_.begin()
        ),
        GAMGInterfaceAgglomerateCoeffs
        (
            fineCoeffs.data(),
            faceRestrictSortAddressing_.data()
        )
    );

    return tcoarseCoeffs;
}


void Foam::GAMGInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE << faceRestrictAddressing_;
}


// ************************************************************************* //
