#include "DeviceConfig.H"

namespace Foam {

    int deviceCount()
    {
        int num_devices;
        CUDA_CALL(cudaGetDeviceCount(&num_devices));
        return num_devices;
    }

    int currentDevice()
    {
        int device;
        CUDA_CALL(cudaGetDevice(&device));
        return device;
    }

    void setCurrentDevice(int device)
    {
        CUDA_CALL(cudaSetDevice(device));
    }

    int deviceComputeCapability(int device)
    {
        cudaDeviceProp deviceProp;
        CUDA_CALL(cudaGetDeviceProperties(&deviceProp, device));
        return 10*deviceProp.major + deviceProp.minor;
    }

    int currentComputeCapability()
    {
        return deviceComputeCapability(currentDevice());
    }

    bool needTextureBind()
    {
        static bool needBind = currentComputeCapability() < 35;
        return needBind;
    }

}
