#pragma once

#include <type_traits>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/unique.h>
#include <thrust/tuple.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>


namespace gpu_api = thrust;

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


namespace Foam
{
    template<>
    inline bool hasAtomicAdd<int>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<long>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<float>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<double>()
    {
        return currentComputeCapability() >= 60;
    }


    template<>
    struct Foam::is_number<int>: std::true_type {};


    template<>
    struct Foam::is_number<long>: std::true_type {};


    template<>
    struct Foam::is_number<float>: std::true_type {};


    template<>
    struct Foam::is_number<double>: std::true_type {};


    template<class T>
    inline bool hasAtomicAdd()
    {
        static_assert(is_number<T>::value, "Number is required");
        return false;
    }
}

