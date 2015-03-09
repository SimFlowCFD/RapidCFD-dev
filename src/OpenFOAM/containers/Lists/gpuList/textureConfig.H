#pragma once

#ifdef CUSP_USE_TEXTURE_MEMORY
    #include <cusp/detail/device/texture.h>
#endif

#ifndef TEXTURE_MINIMUM_SIZE
    #define TEXTURE_MINIMUM_SIZE 200000
#endif

namespace Foam
{

    inline void bind(const scalar* x)
    {
        #ifdef CUSP_USE_TEXTURE_MEMORY
        bind_x(x);
        #endif 
    }

    inline void unbind(const scalar* x)
    {                                    
        #ifdef CUSP_USE_TEXTURE_MEMORY
        unbind_x(x);
        #endif 
    }

    template<bool useTexture>
#ifdef __CUDA_ARCH__ 
    __device__ 
#endif
    inline scalar fetch(const label& i, const scalar* x)
    {
    #ifdef CUSP_USE_TEXTURE_MEMORY
        #ifdef __CUDA_ARCH__
        if(useTexture)
            return fetch_x<true>(i,x);
        else      
        #endif
            return x[i];
    #else
        return x[i];
    #endif
    }
}
