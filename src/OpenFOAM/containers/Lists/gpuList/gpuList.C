#include "gpuList.H"

#include "DeviceMemory.H"
#include "error.H"
#include "pTraits.H"

#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/uninitialized_fill.h>

template<class T>
Foam::gpuList<T>::gpuList(label size)
:
    gpuList(size, T())
{}


template<class T>
Foam::gpuList<T>::gpuList(label size, const T& t)
:
    v_(allocDevice<T>(size)),
    size_(size),
    owner_(true)
{
    thrust::uninitialized_fill(begin(), end(), t);
}


template<class T>
Foam::gpuList<T>::gpuList(const gpuList<T>& list)
:
    gpuList(list.size())
{
    operator=(list);
}


template<class T>
Foam::gpuList<T>::gpuList
(
    const gpuList<T>& list,
    const label subSize
)
:
    v_(list.v_),
    size_(subSize),
    owner_(false)
{}


template<class T>
Foam::gpuList<T>::gpuList
(
    const gpuList<T>& list,
    const label subSize,
    const label startIndex
)
:
    v_(list.v_+startIndex),
    size_(subSize),
    owner_(false)
{}


template<class T>
Foam::gpuList<T>::gpuList(const Xfer<gpuList<T> >& lst)
{
    transfer(lst());
}


template<class T>
Foam::gpuList<T>::gpuList(const UList<T>& list)
:
    gpuList(list.size())
{
    operator=(list);
}


// Construct as copy or re-use as specified.
template<class T>
Foam::gpuList<T>::gpuList(gpuList<T>& a, bool reUse) 
{
    if (reUse)
    {
        v_ = a.v_;
        size_ = a.size_;
        owner_ = a.owner_;

        a.v_ = 0;
        a.size_ = 0;
        a.owner_ = false;
    }
    else
    { 
        v_ = allocDevice<T>(a.size());
        size_ = a.size();
        owner_ = true;
        thrust::uninitialized_fill(begin(), end(), T());
        operator=(a);
    }
}


template<class T>
Foam::gpuList<T>::gpuList()
:
    v_(0),
    size_(0),
    owner_(false)
{}


template<class T>
const Foam::gpuList<T>& Foam::gpuList<T>::null()
{
    return NullObjectRef<gpuList<T>>();
}


template<class T>
Foam::Xfer<Foam::gpuList<T> > Foam::gpuList<T>::xfer()
{
    return xferMove(*this);
}


template<class T>
void Foam::gpuList<T>::transfer(gpuList<T>& a)
{ 
    if (v_ && owner_) freeDevice(v_);

    v_ = a.v_;
    size_ = a.size_;
    owner_ = a.owner_;

    a.v_ = 0;
    a.size_ = 0;
    a.owner_ = false;
}


template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a, label size, label start)
{ 
    if (v_ && owner_) freeDevice(v_);
    v_ = a.v_ + start;
    size_ = size;
    owner_ = false;
}


template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a, label size)
{ 
    setDelegate(a, size, 0);
}


template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a)
{ 
    setDelegate(a, a.size(), 0);
}


template<class T>
Foam::gpuList<T>::~gpuList()
{
    if (v_ && owner_)
    { 
        try
        {
             freeDevice(v_);
        }
        catch(std::runtime_error &e)
        {
            #ifdef FULLDEBUG
                WarningIn("Foam::gpuList<T>::~gpuList()")
                        <<"Error during vector destruction: " << e.what() << endl;
            #endif
        }
    }
}


template<class T>
T Foam::gpuList<T>::first() const
{
    return get(0);
}


template<class T>
T Foam::gpuList<T>::get(const label n) const
{
    if(v_ && n < size_)
    {
        return *const_iterator(v_+n);
    }
    else
    {
        //maybe exception??
        return pTraits<T>::zero;
    }
}


template<class T>
void Foam::gpuList<T>::set(const label n, const T val)
{
    if(v_)
    {
        *iterator(v_+n) = val;
    }
}


template<class T>
void Foam::gpuList<T>::clear()
{
    if(v_ && owner_)
    {
         freeDevice(v_);
    }
    v_ = 0;
    size_ = 0;
    owner_ = false;
}


template<class T>
void Foam::gpuList<T>::setSize(label newSize, const T val)
{
    label oldSize = size_;
    setSize(newSize);
    if(newSize > oldSize)
        thrust::fill(begin()+oldSize, end(), val);
}


template<class T>
void Foam::gpuList<T>::setSize(label newSize)
{
    if (newSize != size_)
    {
        if (newSize > 0)
        {
            T* nv = allocDevice<T>(newSize);
            label span = min(size_, newSize);

            thrust::copy(begin(), begin()+span, iterator(nv));
            thrust::uninitialized_fill(begin()+span, end(), T());

            if (v_ && owner_) freeDevice(v_);

            size_ = newSize;
            v_ = nv;
            owner_ = true;
        }
        else
        {
            clear();
        }
    }
}

template<class T>
void Foam::gpuList<T>::operator=(const T& t)
{
    thrust::fill(begin(), end(), t);
}


template<class T>
void Foam::gpuList<T>::operator=(const gpuList<T>& l)
{
    setSize(l.size());
    thrust::copy(l.begin(), l.end(), begin());
}


template<class T>
void Foam::gpuList<T>::operator=(const UList<T>& l)
{
    setSize(l.size());
    thrust::copy(l.begin(), l.end(), begin());
}


template<class T>
void Foam::gpuList<T>::operator=(const List<T>& l)
{
    operator=(static_cast<const UList<T>&>(l));
}


#include "gpuListIO.C"

