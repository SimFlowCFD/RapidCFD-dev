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

#include "processorGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        processorGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        processorGAMGInterfaceField,
        lduInterfaceField
    );

    struct processorGAMGInterfaceFieldFunctor
    {
         __HOST____DEVICE__
         scalar operator()(const scalar& f,const thrust::tuple<scalar,scalar>& t)
         {
             return f - thrust::get<0>(t)*thrust::get<1>(t);
         }
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorGAMGInterfaceField::processorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const processorLduInterfaceField& p =
        refCast<const processorLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::processorGAMGInterfaceField::processorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorGAMGInterfaceField::~processorGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorGAMGInterfaceField::initInterfaceMatrixUpdate
(
    scalargpuField&,
    const scalargpuField& psiInternal,
    const scalargpuField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = comm();

    scalargpuField tmpf;
    procInterface_.interfaceInternalField(psiInternal, tmpf);

    if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
    {
        // Fast path.
        scalarSendBuf_.setSize(tmpf.size());
        thrust::copy(tmpf.begin(),tmpf.end(),scalarSendBuf_.begin());
        scalarReceiveBuf_.setSize(scalarSendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );
    }
    else
    {
        procInterface_.compressedSend(commsType, tmpf);
    }

    const_cast<processorGAMGInterfaceField&>(*this).updatedMatrix() = false;

    UPstream::warnComm = oldWarn;
}


void Foam::processorGAMGInterfaceField::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField&,
    const scalargpuField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (updatedMatrix())
    {
        return;
    }

    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = comm();

    const labelgpuList& faceCells = procInterface_.faceCells();

    if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;

        // Consume straight from scalarReceiveBuf_
 
        // Transform according to the transformation tensor
        scalargpuField pnf(scalarReceiveBuf_);
        transformCoupleField(pnf, cmpt);

        // Multiply the field by coefficients and add into the result        
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.begin()
            ),
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                coeffs.begin(),
                pnf.begin()
            )),
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.begin()
            ),
            processorGAMGInterfaceFieldFunctor()
        );
    }
    else
    {
        scalarField pnfTmp
        (
            procInterface_.compressedReceive<scalar>(commsType, coeffs.size())
        );
        scalargpuField pnf(pnfTmp);

        transformCoupleField(pnf, cmpt);

        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.begin()
            ),
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                coeffs.begin(),
                pnf.begin()
            )),
            thrust::make_permutation_iterator
            (
                result.begin(),
                faceCells.begin()
            ),
            processorGAMGInterfaceFieldFunctor()
        );
    }

    const_cast<processorGAMGInterfaceField&>(*this).updatedMatrix() = true;

    UPstream::warnComm = oldWarn;
}


// ************************************************************************* //
