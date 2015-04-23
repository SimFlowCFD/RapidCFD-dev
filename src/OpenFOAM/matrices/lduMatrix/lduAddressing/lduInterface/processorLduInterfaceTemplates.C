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

#include "processorLduInterface.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions * * *  * * * * * * * * * * //

namespace Foam
{

template<class To, class From, bool compress>
struct compressFunctor
{
    To* to;
    const From* from;
    const scalar* slast;
    const label nCmpts;

    compressFunctor
    (
        To* _to,
        const From* _from,
        const scalar* _slast,
        label _nCmpts
    ):
        to(_to),
        from(_from),
        slast(_slast),
        nCmpts(_nCmpts)
    {}

    __HOST____DEVICE__
    void operator()(const label& i)
    {
        if(compress)
            to[i] = (To) (from[i] - slast[i%nCmpts]);
        else
            to[i] = (To) (from[i] + slast[i%nCmpts]);
    }
};

}

template<class Type>
void Foam::processorLduInterface::send
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    label nBytes = f.byteSize();

    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        OPstream::write
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<const char*>(f.begin()),
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        resizeBuf(receiveBuf_, nBytes);

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            receiveBuf_.begin(),
            nBytes,
            tag(),
            comm()
        );

        resizeBuf(sendBuf_, nBytes);
        memcpy(sendBuf_.begin(), f.begin(), nBytes);

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            sendBuf_.begin(),
            nBytes,
            tag(),
            comm()
        );
    }
    else
    {
        FatalErrorIn("processorLduInterface::send")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}

template<class Type>
void Foam::processorLduInterface::send
(
    const Pstream::commsTypes commsType,
    const gpuList<Type>& f
) const
{
    label nBytes = f.byteSize();

    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        const char* sendData;
        if(Pstream::gpuDirectTransfer)
        {
            sendData = reinterpret_cast<const char*>(f.data());
        }
        else
        {
            resizeBuf(sendBuf_, nBytes);
            CUDA_CALL(cudaMemcpy(sendBuf_.begin(), f.data(), nBytes,cudaMemcpyDeviceToHost));
            sendData = sendBuf_.begin();
        }

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            sendData,
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        char* receive;
        const char* send;

        if(Pstream::gpuDirectTransfer)
        {
            resizeBuf(gpuReceiveBuf_, nBytes);
            resizeBuf(gpuSendBuf_, nBytes);

            CUDA_CALL(cudaMemcpy(gpuSendBuf_.data(), f.data(), nBytes,cudaMemcpyDeviceToDevice));

            send = gpuSendBuf_.data();
            receive = gpuReceiveBuf_.data();
        }
        else
        {
            resizeBuf(receiveBuf_, nBytes);
            resizeBuf(sendBuf_, nBytes);

            CUDA_CALL(cudaMemcpy(sendBuf_.begin(), f.data(), nBytes,cudaMemcpyDeviceToHost));

            send = sendBuf_.begin();
            receive = receiveBuf_.begin();
        }

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            receive,
            nBytes,
            tag(),
            comm()
        );

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            send,
            nBytes,
            tag(),
            comm()
        );
    }
    else
    {
        FatalErrorIn("processorLduInterface::send")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}


template<class Type>
void Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        IPstream::read
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<char*>(f.begin()),
            f.byteSize(),
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        memcpy(f.begin(), receiveBuf_.begin(), f.byteSize());
    }
    else
    {
        FatalErrorIn("processorLduInterface::receive")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}

template<class Type>
void Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    gpuList<Type>& f
) const
{
    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        char * read;
        if(Pstream::gpuDirectTransfer)
        {
            read = reinterpret_cast<char*>(f.data());
        }
        else
        {
            resizeBuf(receiveBuf_, f.byteSize());
            read = receiveBuf_.begin();
        }

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<char*>(read),
            f.byteSize(),
            tag(),
            comm()
        );

        if( ! Pstream::gpuDirectTransfer)
        {
            CUDA_CALL(cudaMemcpy(f.data(), receiveBuf_.data(), f.byteSize(),cudaMemcpyHostToDevice));
        }
    }
    else if (commsType == Pstream::nonBlocking)
    {
        if(Pstream::gpuDirectTransfer)
        {
            CUDA_CALL(cudaMemcpy(f.data(), gpuReceiveBuf_.data(), f.byteSize(),cudaMemcpyDeviceToDevice));
        }
        else
        {
            CUDA_CALL(cudaMemcpy(f.data(), receiveBuf_.data(), f.byteSize(),cudaMemcpyHostToDevice));
        }
    }
    else
    {
        FatalErrorIn("processorLduInterface::receive")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<gpuField<Type> > tf(new gpuField<Type>(size));
    receive(commsType, tf());
    return tf;
}

template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];
        resizeBuf(sendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(sendBuf_.begin());

        for (register label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        reinterpret_cast<Type&>(fArray[nm1]) = f.last();

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == Pstream::nonBlocking)
        {
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else
        {
            FatalErrorIn("processorLduInterface::compressedSend")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }
    }
    else
    {
        this->send(commsType, f);
    }
}

template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const Pstream::commsTypes commsType,
    const gpuList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.data());
        const scalar *slast = &sArray[nm1];
        resizeBuf(gpuSendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(gpuSendBuf_.data());

        thrust::for_each
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(nm1),
            compressFunctor<float,scalar,true>
            (
                fArray,
                sArray,
                slast,
                nCmpts
             )
        );

        CUDA_CALL(cudaMemcpy(fArray+nm1, f.data() + (f.size() - 1), sizeof(Type), cudaMemcpyDeviceToDevice));

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            const char* sendData;
            if(Pstream::gpuDirectTransfer)
            {
                sendData = gpuSendBuf_.data();
            }
            else
            {
                resizeBuf(sendBuf_, nBytes);
                CUDA_CALL(cudaMemcpy(sendBuf_.begin(), gpuSendBuf_.data(), nBytes,cudaMemcpyDeviceToHost));
                sendData = sendBuf_.begin();
            }

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendData,
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == Pstream::nonBlocking)
        {
            const char* sendData;
            char * readData;

            if(Pstream::gpuDirectTransfer)
            {
                resizeBuf(gpuReceiveBuf_, nBytes);
                readData = gpuReceiveBuf_.data();
                sendData = gpuSendBuf_.data();
            }
            else
            {
                resizeBuf(receiveBuf_, nBytes);
                resizeBuf(sendBuf_, nBytes);

                CUDA_CALL(cudaMemcpy(sendBuf_.begin(), gpuSendBuf_.data(), nBytes,cudaMemcpyDeviceToHost));

                sendData = sendBuf_.begin();
                readData = receiveBuf_.begin();
            }

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                readData,
                nBytes,
                tag(),
                comm()
            );

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendData,
                nBytes,
                tag(),
                comm()
            );
        }
        else
        {
            FatalErrorIn("processorLduInterface::compressedSend")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }
    }
    else
    {
        this->send(commsType, f);
    }
}

template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType != Pstream::nonBlocking)
        {
            FatalErrorIn("processorLduInterface::compressedReceive")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }

        const float *fArray =
            reinterpret_cast<const float*>(receiveBuf_.begin());
        f.last() = reinterpret_cast<const Type&>(fArray[nm1]);
        scalar *sArray = reinterpret_cast<scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];

        for (register label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }
    }
    else
    {
        this->receive<Type>(commsType, f);
    }
}

template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    gpuList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        resizeBuf(gpuReceiveBuf_, nBytes);
        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            char* readData;
            if(Pstream::gpuDirectTransfer)
            {
                readData = gpuReceiveBuf_.data();
            }
            else
            {
                resizeBuf(receiveBuf_, nBytes);
                readData = receiveBuf_.begin();
            }

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                readData,
                nBytes,
                tag(),
                comm()
            );

            if( ! Pstream::gpuDirectTransfer)
            {
                CUDA_CALL(cudaMemcpy(gpuReceiveBuf_.data(), receiveBuf_.data(), nBytes,cudaMemcpyHostToDevice));
            }
        }
        else if (commsType == Pstream::nonBlocking)
        {
            if( ! Pstream::gpuDirectTransfer)
            {
                CUDA_CALL(cudaMemcpy(gpuReceiveBuf_.data(), receiveBuf_.data(), nBytes,cudaMemcpyHostToDevice));
            }
        }
        else
        {
            FatalErrorIn("processorLduInterface::compressedReceive")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }

        const float *fArray =
            reinterpret_cast<const float*>(gpuReceiveBuf_.data());

        CUDA_CALL(cudaMemcpy(f.data()+(f.size() - 1),fArray+nm1, sizeof(Type), cudaMemcpyDeviceToDevice));

        scalar *sArray = reinterpret_cast<scalar*>(f.data());
        const scalar *slast = &sArray[nm1];

        thrust::for_each
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(nm1),
            compressFunctor<scalar,float,false>
            (
                sArray,
                fArray,
                slast,
                nCmpts
             )
        );
    }
    else
    {
        this->receive<Type>(commsType, f);
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<gpuField<Type> > tf(new gpuField<Type>(size));
    compressedReceive(commsType, tf());
    return tf;
}


// ************************************************************************* //
