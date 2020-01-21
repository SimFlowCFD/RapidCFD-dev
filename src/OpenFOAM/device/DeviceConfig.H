#pragma once

#include <stdio.h>
#include <type_traits>

#define CUDA_CALL(x) do { if((x) != cudaSuccess) {         \
 printf("Error at %s:%d\n",__FILE__,__LINE__);             \
 printf("%s\n",cudaGetErrorString(cudaPeekAtLastError())); \
 cudaDeviceReset();                                        \
 ::exit(static_cast<int>(cudaPeekAtLastError()));          \
 }} while(0)

#define GPU_ERROR_CHECK()                                  \
 cudaDeviceSynchronize();                                  \
 CUDA_CALL( cudaPeekAtLastError());

#define GPU_ERROR_CHECK_ASYNC()                            \
 CUDA_CALL(cudaPeekAtLastError());

namespace Foam
{
    int deviceCount();
    int currentDevice();
    void setCurrentDevice(int device);
    int deviceComputeCapability(int device);
    int currentComputeCapability();
    bool needTextureBind();
    
    template<class T>
    bool hasAtomicAdd();

    template<class T>
    struct is_number: std::false_type {};
}

#include "DeviceConfigI.H"

