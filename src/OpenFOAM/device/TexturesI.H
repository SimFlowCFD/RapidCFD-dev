#pragma once

namespace Foam
{

#if __CUDA_ARCH__ >= 350

template<class T>
inline __device__ T textures<T>::operator[](const int& i) const
{
    return __ldg(data + i);
}

#else

template<>
inline __device__ int textures<int>::operator[](const int& i) const
{
    return tex1Dfetch<int>(tex, i);
}

template<>
inline __device__ long textures<long>::operator[](const int& i) const
{
    int2 v = tex1Dfetch<int2>(tex, i);
    long res = v.y;
    res = res << 32;
    return res | v.x;
}

template<>
inline __device__ float textures<float>::operator[](const int& i) const
{
    return tex1Dfetch<float>(tex, i);
}

template<>
inline __device__ double textures<double>::operator[](const int& i) const
{
    int2 v = tex1Dfetch<int2>(tex, i);
    return __hiloint2double(v.y, v.x);
}

#endif

template<class T>
inline void textureBind<T>::bind(int n)
{
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(cudaResourceDesc));

    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = const_cast<T*>(data);
    resDesc.res.linear.sizeInBytes = n*sizeof(T);

    initResourceDesc(resDesc);

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(cudaTextureDesc));
    texDesc.readMode = cudaReadModeElementType;

    cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
}

template<>
inline void textureBind<int>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
}

template<>
inline void textureBind<long>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
}

template<>
inline void textureBind<float>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
}

template<>
inline void textureBind<double>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
}

template<class T>
inline void textureBind<T>::initResourceDesc(cudaResourceDesc& resDesc)
{
    static_assert(is_number<T>::value, "Only numbers can be used as texture resources");
}


}
