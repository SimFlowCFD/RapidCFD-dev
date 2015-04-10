#pragma once

#include "gpuList.H"

namespace Foam
{

template<class T>
struct textures
{
private:
    cudaTextureObject_t tex;

    inline void initResourceDesc(cudaResourceDesc& resDesc);
    void init(int n, T* data);

public:
    textures(int n, T* data):
        tex(0)
    {
        init(n,data);
    }

    textures(const gpuList<T>& list):
        tex(0)
    {
        init(list.size(),const_cast<T*>(list.data()));
    }

    inline __device__ T operator[](const int& i);

    void destroy()
    {
        cudaDestroyTextureObject(tex);
    }
};

template<class T>
inline void textures<T>::initResourceDesc(cudaResourceDesc& resDesc)
{
}

template<class T>
inline void textures<T>::init(int n, T* data)
{
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(cudaResourceDesc));
   
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = data;
    resDesc.res.linear.sizeInBytes = n*sizeof(T);

    initResourceDesc(resDesc); 

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(cudaTextureDesc));
    texDesc.readMode = cudaReadModeElementType;

    cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
}

template<>
inline void textures<float>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32; 
}


template<>
inline void textures<int>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32; 
}

template<>
inline void textures<double>::initResourceDesc(cudaResourceDesc& resDesc)
{
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32; 
    resDesc.res.linear.desc.y = 32; 
}

template<>
inline __device__ float textures<float>::operator[](const int& i)
{
    return tex1Dfetch<float>(tex, i);
}

template<>
inline __device__ int textures<int>::operator[](const int& i)
{
    return tex1Dfetch<int>(tex, i);
}

template<>
inline __device__ double textures<double>::operator[](const int& i)
{
    int2 v = tex1Dfetch<int2>(tex, i);
    return __hiloint2double(v.y, v.x);
}

}
