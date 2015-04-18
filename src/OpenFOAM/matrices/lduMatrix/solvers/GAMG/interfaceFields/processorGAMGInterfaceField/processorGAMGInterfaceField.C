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
#include "GAMGInterfaceFunctors.H"

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

    procInterface_.interfaceInternalField(psiInternal, scalargpuSendBuf_);

    if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
    {
        std::streamsize nBytes = scalargpuSendBuf_.byteSize();

        scalar* readData;
        const scalar* sendData;

        if(Pstream::gpuDirectTransfer)
        {
            // Fast path.
            scalargpuReceiveBuf_.setSize(scalargpuSendBuf_.size());

            sendData = scalargpuSendBuf_.data();
            readData = scalargpuReceiveBuf_.data();
        }
        else
        {
            scalarSendBuf_.setSize(scalargpuSendBuf_.size());
            scalarReceiveBuf_.setSize(scalarSendBuf_.size());
            thrust::copy
            (
                scalargpuSendBuf_.begin(),
                scalargpuSendBuf_.end(),
                scalarSendBuf_.begin()
            );

            sendData = scalarSendBuf_.begin();
            readData = scalarReceiveBuf_.begin();
        }

        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(readData),
            nBytes,
            procInterface_.tag(),
            comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
        (
            Pstream::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(sendData),
            nBytes,
            procInterface_.tag(),
            comm()
        );
    }
    else
    {
        procInterface_.compressedSend(commsType, scalargpuSendBuf_);
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

        if( ! Pstream::gpuDirectTransfer)
        {
            scalargpuReceiveBuf_ = scalarReceiveBuf_;
        }

        // Transform according to the transformation tensor
        transformCoupleField(scalargpuReceiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result  
        GAMGUpdateInterfaceMatrix
        (
            result,
            coeffs,
            scalargpuReceiveBuf_,
            procInterface_
        );      
    }
    else
    {
        scalargpuReceiveBuf_.setSize(coeffs.size());
        procInterface_.compressedReceive<scalar>(commsType, scalargpuReceiveBuf_);

        transformCoupleField(scalargpuReceiveBuf_, cmpt);

        GAMGUpdateInterfaceMatrix
        (
            result,
            coeffs,
            scalargpuReceiveBuf_,
            procInterface_
        ); 
    }

    const_cast<processorGAMGInterfaceField&>(*this).updatedMatrix() = true;

    UPstream::warnComm = oldWarn;
}


// ************************************************************************* //
