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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p)),
    gpuSendBuf_(0),
    gpuReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(0),
    scalargpuReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p)),
    gpuSendBuf_(0),
    gpuReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(0),
    scalargpuReceiveBuf_(0)
{}


// Construct by mapping given processorFvPatchField<Type>
template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p)),
    gpuSendBuf_(0),
    gpuReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(0),
    scalargpuReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
    if (debug && !ptf.ready())
    {
        FatalErrorIn("processorFvPatchField<Type>::processorFvPatchField(..)")
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p)),
    gpuSendBuf_(0),
    gpuReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(0),
    scalargpuReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    gpuSendBuf_(ptf.gpuSendBuf_.xfer()),
    gpuReceiveBuf_(ptf.gpuReceiveBuf_.xfer()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(ptf.scalargpuSendBuf_.xfer()),
    scalargpuReceiveBuf_(ptf.scalargpuReceiveBuf_.xfer())
{
    if (debug && !ptf.ready())
    {
        FatalErrorIn("processorFvPatchField<Type>::processorFvPatchField(..)")
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    gpuSendBuf_(0),
    gpuReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalargpuSendBuf_(0),
    scalargpuReceiveBuf_(0)
{
    if (debug && !ptf.ready())
    {
        FatalErrorIn("processorFvPatchField<Type>::processorFvPatchField(..)")
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::~processorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::processorFvPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorIn("processorFvPatchField<Type>::patchNeighbourField()")
            << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        this->patchInternalField(gpuSendBuf_);

        if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
        {
            std::streamsize nBytes = gpuSendBuf_.byteSize();

            Type* receive;
            const Type* send;

            this->setSize(gpuSendBuf_.size());
            if(Pstream::gpuDirectTransfer)
            {
                // Fast path.
                send = gpuSendBuf_.data();
                receive = this->data();
            }
            else
            {
                sendBuf_.setSize(gpuSendBuf_.size());
                receiveBuf_.setSize(sendBuf_.size());
                thrust::copy
                (
                    gpuSendBuf_.begin(),
                    gpuSendBuf_.end(),
                    sendBuf_.begin()
                );

                send = sendBuf_.begin();
                receive = receiveBuf_.begin();
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
            procPatch_.compressedSend(commsType, gpuSendBuf_);
        }
    }
}


template<class Type>
void Foam::processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
        {
            // Fast path. Received into *this

            if
            (
                outstandingRecvRequest_ >= 0
             && outstandingRecvRequest_ < Pstream::nRequests()
            )
            {
                UPstream::waitRequest(outstandingRecvRequest_);
            }
            outstandingSendRequest_ = -1;
            outstandingRecvRequest_ = -1;

            if( ! Pstream::gpuDirectTransfer)
            {
                scalargpuReceiveBuf_ = scalarReceiveBuf_;
                thrust::copy
                (
                    receiveBuf_.begin(),
                    receiveBuf_.end(),
                    this->begin()
                );
            }
        }
        else
        {
            procPatch_.compressedReceive<Type>(commsType, *this);
        }
       
        if (doTransform())
        {
            transform(*this, procPatch_.getForwardT(), *this);
        }
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::processorFvPatchField<Type>::snGrad
(
    const scalargpuField& deltaCoeffs
) const
{
    return deltaCoeffs*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
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
                "processorFvPatchField<Type>::initInterfaceMatrixUpdate(..)"
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

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField&,
    const scalargpuField& coeffs,
    const direction cmpt,
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

        // Consume straight from scalarReceiveBuf_

        if( ! Pstream::gpuDirectTransfer)
        {
            scalargpuReceiveBuf_ = scalarReceiveBuf_;
        }

        // Transform according to the transformation tensor
        transformCoupleField(scalargpuReceiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result
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

        // Transform according to the transformation tensor
        transformCoupleField(scalargpuReceiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result
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

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    gpuField<Type>&,
    const gpuField<Type>& psiInternal,
    const scalargpuField&,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, gpuSendBuf_);

    if (commsType == Pstream::nonBlocking && !Pstream::floatTransfer)
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorIn
            (
                "processorFvPatchField<Type>::initInterfaceMatrixUpdate(..)"
            )   << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        std::streamsize nBytes = gpuSendBuf_.byteSize();

        Type* receive;
        const Type* send;

        if(Pstream::gpuDirectTransfer)
        {
            // Fast path.
            gpuReceiveBuf_.setSize(gpuSendBuf_.size());

            send = gpuSendBuf_.data();
            receive = gpuReceiveBuf_.data();
        }
        else
        {
            sendBuf_.setSize(gpuSendBuf_.size());
            receiveBuf_.setSize(sendBuf_.size());
            thrust::copy
            (
                gpuSendBuf_.begin(),
                gpuSendBuf_.end(),
                sendBuf_.begin()
            );

            send = sendBuf_.begin();
            receive = receiveBuf_.begin();
        }

        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(receive),
            nBytes,
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
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
        procPatch_.compressedSend(commsType, gpuSendBuf_);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    gpuField<Type>& result,
    const gpuField<Type>&,
    const scalargpuField& coeffs,
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

        // Consume straight from receiveBuf_

        if( ! Pstream::gpuDirectTransfer)
        {
            gpuReceiveBuf_ = receiveBuf_;
        }

        // Transform according to the transformation tensor
        transformCoupleField(gpuReceiveBuf_);

        // Multiply the field by coefficients and add into the result
        matrixPatchOperation
        (
            this->patch().index(),
            result,
            this->patch().boundaryMesh().mesh().lduAddr(),
            matrixInterfaceFunctor<Type>
            (
                coeffs.data(),
                gpuReceiveBuf_.data()
            )
        );
    }
    else
    {
        gpuReceiveBuf_.setSize(this->size());
        procPatch_.compressedReceive<Type>(commsType, gpuReceiveBuf_);

        // Transform according to the transformation tensor
        transformCoupleField(gpuReceiveBuf_);

        // Multiply the field by coefficients and add into the result
        matrixPatchOperation
        (
            this->patch().index(),
            result,
            this->patch().boundaryMesh().mesh().lduAddr(),
            matrixInterfaceFunctor<Type>
            (
                coeffs.data(),
                gpuReceiveBuf_.data()
            )
        );
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
bool Foam::processorFvPatchField<Type>::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingSendRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingRecvRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}


// ************************************************************************* //
