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
        OPstream::write
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<const char*>(f.data()),
            nBytes,
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        resizeBuf(gpuReceiveBuf_, nBytes);

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            gpuReceiveBuf_.data(),
            nBytes,
            tag(),
            comm()
        );

        resizeBuf(gpuSendBuf_, nBytes);
        cudaMemcpy(gpuSendBuf_.data(), f.data(), nBytes,cudaMemcpyDeviceToDevice);

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            gpuSendBuf_.data(),
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
        IPstream::read
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<char*>(f.data()),
            f.byteSize(),
            tag(),
            comm()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        cudaMemcpy(f.data(), gpuReceiveBuf_.data(), f.byteSize(),cudaMemcpyDeviceToDevice);
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

        cudaMemcpy(fArray+nm1, f.data() + (f.size() - 1), sizeof(Type), cudaMemcpyDeviceToDevice);

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            OPstream::write
            (
                commsType,
                neighbProcNo(),
                gpuSendBuf_.data(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == Pstream::nonBlocking)
        {
            resizeBuf(gpuReceiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                gpuReceiveBuf_.data(),
                nBytes,
                tag(),
                comm()
            );

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                gpuSendBuf_.data(),
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

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            resizeBuf(gpuReceiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                gpuReceiveBuf_.data(),
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
            reinterpret_cast<const float*>(gpuReceiveBuf_.data());

        cudaMemcpy(f.data()+(f.size() - 1),fArray+nm1, sizeof(Type), cudaMemcpyDeviceToDevice);

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
