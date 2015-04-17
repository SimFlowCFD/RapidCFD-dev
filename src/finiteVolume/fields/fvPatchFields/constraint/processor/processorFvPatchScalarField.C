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

#include "processorFvPatchScalarField.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    scalargpuField&,
    const scalargpuField& psiInternal,
    const scalargpuField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, scalargpuSendBuf_);

    if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorIn
            (
                "processorFvPatchField<scalar>::initInterfaceMatrixUpdate(..)"
            )   << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        std::streamsize nBytes = scalargpuSendBuf_.byteSize();

        scalar* receive;
        const scalar* send;

        if(Pstream::gpuDirectTransfer)
        {
            // Fast path.
            scalargpuReceiveBuf_.setSize(scalargpuSendBuf_.size());

            send = scalargpuSendBuf_.data();
            receive = scalargpuReceiveBuf_.data();
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

            send = scalarSendBuf_.begin();
            receive = scalarReceiveBuf_.begin();
        }

        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(receive),
            nBytes,
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(send),
            nBytes,
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, scalargpuSendBuf_);
    }

    const_cast<processorFvPatchField<scalar>&>(*this).updatedMatrix() = false;
}


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField&,
    const scalargpuField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

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

        if( ! Pstream::gpuDirectTransfer)
        {
            scalargpuReceiveBuf_ = scalarReceiveBuf_;
        }

        // Consume straight from scalargpuReceiveBuf_
        matrixPatchOperation
        (
            this->patch().index(),
            result,
            this->patch().boundaryMesh().mesh().lduAddr(),
            matrixInterfaceFunctor<scalar>
            (
                coeffs.data(),
                scalargpuReceiveBuf_.data()
            )
       );
    }
    else
    {
        scalargpuReceiveBuf_.setSize(this->size());
        procPatch_.compressedReceive<scalar>(commsType, scalargpuReceiveBuf_);

        matrixPatchOperation
        (
            this->patch().index(),
            result,
            this->patch().boundaryMesh().mesh().lduAddr(),
            matrixInterfaceFunctor<scalar>
            (
                coeffs.data(),
                scalargpuReceiveBuf_.data()
            )
        );
    }

    const_cast<processorFvPatchField<scalar>&>(*this).updatedMatrix() = true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
