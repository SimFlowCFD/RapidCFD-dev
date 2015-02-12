#ifndef gpuConfig_H
#define gpuConfig_H

#if WM_GPU == CUDA

#define GPU_FUNCTOR(T) T
#define GPU_TEMPLATE_FUNCTOR(T) T

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


#define gpuErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define GPU_ERROR_CHECK()                        \
cudaDeviceSynchronize();                         \
gpuErrorCheck( cudaPeekAtLastError() )    

namespace Foam
{

inline void gpuAssert(cudaError_t code, const char *file, int line)
{
   if (code != cudaSuccess) 
   {
      
      Info << "GPUassert: " << cudaGetErrorString(code)
           << ", file: " << file
           << ", line: " << line << endl;

      cudaDeviceReset();
      ::exit(static_cast<label>(code));
   }
}

inline int getGpuDeviceCount()
{
    int num_devices;
    cudaGetDeviceCount(&num_devices);
    return num_devices;
}

inline void setGpuDevice(int device)
{
   cudaSetDevice(device);
}

}

#else
#error "Currently only CUDA is supported."
#endif

#endif
