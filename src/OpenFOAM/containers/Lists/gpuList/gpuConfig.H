#ifndef gpuConfig_H
#define gpuConfig_H

#ifdef __CUDACC__

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>

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

inline int getGpuDeviceCount()
{
    int num_devices;
    CUDA_CALL(cudaGetDeviceCount(&num_devices));
    return num_devices;
}

inline void setGpuDevice(int device)
{
   CUDA_CALL(cudaSetDevice(device));
}

}

#else
#error "Currently only CUDA is supported."
#endif

#endif
