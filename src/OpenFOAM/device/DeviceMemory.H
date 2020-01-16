#pragma once

#include "DeviceConfig.H"

namespace Foam {

    template<class T>
    inline void copyHostToDevice(T* dst, const T* src, const label size)
    {
        CUDA_CALL(cudaMemcpy(dst, src, size*sizeof(T), cudaMemcpyHostToDevice));
    }

    template<class T>
    inline void copyDeviceToHost(T* dst, const T* src, const label size)
    {
        CUDA_CALL(cudaMemcpy(dst, src, size*sizeof(T), cudaMemcpyDeviceToHost));
    }

    template<class T>
    inline void copyDeviceToDevice(T* dst, const T* src, const label size)
    {
        CUDA_CALL(cudaMemcpy(dst, src, size*sizeof(T), cudaMemcpyDeviceToDevice));
    }

    template<class T>
    inline T* allocPageLocked(const label size)
    {
        T* ptr;
        CUDA_CALL(cudaMallocHost(&ptr, size*sizeof(T)));
        return ptr;
    }

    template<class T>
    inline void freePageLocked(T* ptr)
    {
         CUDA_CALL(cudaFreeHost(ptr));
    }

}



