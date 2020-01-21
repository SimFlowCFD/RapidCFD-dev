#pragma once

namespace Foam {
    template<class T>
    inline T* allocDevice(const label size);
    template<class T>
    inline void freeDevice(T* ptr);

    template<class T>
    inline T* allocPageLocked(const label size);
    template<class T>
    inline void freePageLocked(T* ptr);

    inline void copyHostToDevice(void* dst, const void* src, const label size);
    inline void copyDeviceToHost(void* dst, const void* src, const label size);
    inline void copyDeviceToDevice(void* dst, const void* src, const label size);
}

#include "DeviceMemoryI.H"

