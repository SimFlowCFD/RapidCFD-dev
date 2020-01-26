#pragma once

#include "gpuList.H"
#include "DeviceConfig.H"

namespace Foam
{

template<class T>
class textureBind;

template<class T>
class textures
{
    cudaTextureObject_t tex;
    const T* data;

    textures(cudaTextureObject_t _tex, const T* _data):
        tex(_tex), data(_data) {}

public:

    friend class textureBind<T>;

    inline __device__ T operator[](const int& i) const;
};


template<class T>
class textureBind
{
private:
    cudaTextureObject_t tex;
    const T* data;

    inline void initResourceDesc(cudaResourceDesc& resDesc);
    void bind(int n);

public:

    textureBind(int n, T* _data):
        tex(0),
        data(_data)
    {
        if(needTextureBind())
            bind(n);
    }

    textureBind(const gpuList<T>& list):
        tex(0),
        data(list.data())
    {
        if(needTextureBind())
            bind(list.size());
    }

    ~textureBind()
    {
        if(needTextureBind())
            cudaDestroyTextureObject(tex);
    }

    textures<T> operator()() const
    {
        return textures<T>(tex, data);
    }
};

}

#include "TexturesI.H"

