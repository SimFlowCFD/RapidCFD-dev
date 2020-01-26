

inline Foam::DeviceStream::DeviceStream()
{
    CUDA_CALL(cudaStreamCreate(&stream_));
}


inline Foam::DeviceStream::~DeviceStream()
{
    CUDA_CALL(cudaStreamDestroy(stream_));
}


inline void Foam::DeviceStream::synchronize() const
{
    CUDA_CALL(cudaStreamSynchronize(stream_));
}

inline cudaStream_t Foam::DeviceStream::operator()() const
{
    return stream_;
}
