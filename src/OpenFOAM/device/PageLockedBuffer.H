#pragma once

#include "DeviceMemory.H"
#include "UList.H"
#include "Field.H"

namespace Foam
{
template<class T>
inline void makePageLocked(UList<T>&, const label);
template<class T>
inline void clearPageLocked(UList<T>&);

template<class T>
class PageLockedBuffer
{
     mutable Field<T> buffer_;

     PageLockedBuffer(const PageLockedBuffer&) = delete;

public:
    PageLockedBuffer();
    PageLockedBuffer(const label);
    ~PageLockedBuffer();

     Field<T>& buffer(const label) const;
};

}

#include "PageLockedBufferI.H"
