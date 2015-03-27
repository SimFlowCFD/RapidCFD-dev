#include "error.H"
#include "pTraits.H"


template<class T>
inline Foam::gpuList<T>::gpuList()
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    v_ = new gpu_api::device_vector<T>(0);
}

template<class T>
inline Foam::gpuList<T>::gpuList(label size)
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    v_ = new gpu_api::device_vector<T>(size);
}

template<class T>
inline Foam::gpuList<T>::gpuList(label size, const T& t)
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    v_ = new gpu_api::device_vector<T>(size,t);
}

template<class T>
inline Foam::gpuList<T>::gpuList(const gpuList<T>& list)
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    v_ = new gpu_api::device_vector<T>(list.size());
    gpu_api::copy(list.begin(),list.end(),begin());
}

template<class T>
template<class InputIterator>
inline Foam::gpuList<T>::gpuList(InputIterator first, InputIterator last)
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    v_ = new gpu_api::device_vector<T>(last-first);
    gpu_api::copy(first,last,begin());
}

template<class T>
inline Foam::gpuList<T>::gpuList
(
    const gpuList<T>& list,
    const label subSize
)
:
    v_(0),
    size_(subSize),
    start_(0),
    delegate_(0)
{
    delegate_ = const_cast<gpuList<T>*>(&list);
}

template<class T>
inline Foam::gpuList<T>::gpuList
(
    const gpuList<T>& list,
    const label subSize,
    const label startIndex
)
:
    v_(0),
    size_(subSize),
    start_(startIndex),
    delegate_(0)
{
    delegate_ = const_cast<gpuList<T>*>(&list);
}

template<class T>
inline Foam::gpuList<T>::gpuList(const Xfer<gpuList<T> >& lst)
{
    transfer(lst());
}

template<class T>
inline Foam::gpuList<T>::gpuList(const UList<T>& list)
:
    v_(new gpu_api::device_vector<T>(list.size())),
    size_(0),
    start_(0),
    delegate_(0)
{
    this->operator=(list);
}

// Construct as copy or re-use as specified.
template<class T>
Foam::gpuList<T>::gpuList(gpuList<T>& a, bool reUse)
:
    v_(0),
    size_(0),
    start_(0),
    delegate_(0)
{
    if (reUse)
    {
        this->v_ = a.v_;
        this->size_ = a.size_;
        this->start_ = a.start_;
        this->delegate_ = a.delegate_;
        a.v_ = 0;
        a.size_ = 0;
        a.start_ = 0;
        a.delegate_ = 0;
    }
    else
    { 
        this->v_ = new gpu_api::device_vector<T>(a.size());

        this->operator=(a);
    }
}




template<class T>
inline const Foam::gpuList<T>& Foam::gpuList<T>::null()
{
    return *reinterpret_cast< gpuList<T>* >(0);
}

template<class T>
inline Foam::Xfer<Foam::gpuList<T> > Foam::gpuList<T>::xfer()
{
    return xferMove(*this);
}

template<class T>
void Foam::gpuList<T>::transfer(gpuList<T>& a)
{ 
    if (this->v_) delete this->v_;
    this->v_ = 0;

    if(a.v_)
    {
        this->v_ = a.v_;
        a.v_ = 0;
    }
    else
    {
        this->delegate_ = a.delegate_;
        this->size_ = a.size_;
        this->start_ = a.start_;
    }
}

template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a, label size, label start)
{ 
    if (this->v_) delete this->v_;
    this->v_ = 0;

    this->delegate_ = &a;
    this->size_ = size;
    this->start_ = start;
}

template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a, label size)
{ 
    setDelegate(a,size,0);
}

template<class T>
void Foam::gpuList<T>::setDelegate(gpuList<T>& a)
{ 
    setDelegate(a,a.size(),0);
}

template<class T>
inline T* Foam::gpuList<T>::data()
{
    if(v_)
    {
        return gpu_api::raw_pointer_cast(v_->data());
    }
    else
    {
        return delegate_->data()+start_;
    }
}

template<class T>
inline const T* Foam::gpuList<T>::data() const
{
    if(v_)
    {
         return gpu_api::raw_pointer_cast(v_->data());
    }
    else
    {
         return delegate_->data()+start_;
    }
}

template<class T>
inline const Foam::gpuList<T>::iterator Foam::gpuList<T>::begin()
{
    if(v_)
    {
        return v_->begin();
    }
    else
    { 
        iterator it = delegate_->v_->begin();
        return it+start_;
    }
}

template<class T>
inline const Foam::gpuList<T>::iterator Foam::gpuList<T>::end()
{
    if(v_)
    {
         return v_->end();
    }
    else
    {
         iterator it = delegate_->v_->begin();
         return it+(start_+size_);
    }
}

template<class T>
inline const Foam::gpuList<T>::const_iterator Foam::gpuList<T>::begin() const
{
    if(v_)
    {
         return v_->begin();
    }
    else
    {
         const_iterator it = delegate_->v_->begin();
         return it+start_;
    }
}

template<class T>
inline const Foam::gpuList<T>::const_iterator Foam::gpuList<T>::end() const
{
    if(v_)
    {
         return v_->end();
    }
    else
    {
         const_iterator it = delegate_->v_->begin();
         return it+(start_+size_);
    }
}

template<class T>
inline const Foam::gpuList<T>::reverse_iterator Foam::gpuList<T>::rbegin()
{
    if(v_)
    {
         return v_->rbegin();
    }
    else
    {
         reverse_iterator it = delegate_->v_->rbegin();
         return it+(delegate_->size()-start_-size_);
    }
}

template<class T>
inline const Foam::gpuList<T>::reverse_iterator Foam::gpuList<T>::rend()
{
    if(v_)
    {
         return v_->rend();
    }
    else
    {
         reverse_iterator it = delegate_->v_->rbegin();
         return it+(delegate_->size()-start_);
    }
}

template<class T>
inline const Foam::gpuList<T>::const_reverse_iterator Foam::gpuList<T>::rbegin() const
{
    if(v_)
    {
        return v_->rbegin();
    }
    else
    {
        const_reverse_iterator it = delegate_->v_->rbegin();
        return it+(delegate_->size()-start_-size_);
    }
}

template<class T>
inline const Foam::gpuList<T>::const_reverse_iterator Foam::gpuList<T>::rend() const
{
    if(v_)
    {
        return v_->rend();
    }
    else
    {
        const_reverse_iterator it = delegate_->v_->rbegin();
        return it+(delegate_->size()-start_);
    }
}


template<class T>
inline T Foam::gpuList<T>::first() const
{
    if(v_)
    {
        if(v_->size())
        {
             return (*v_)[0];
        }
        else
        {
            //maybe exception??
            return pTraits<T>::zero;
        }
    }
    else
    {
        return delegate_->get(start_);
    }
}

template<class T>
inline T Foam::gpuList<T>::get(const label n) const
{
    if(v_)
    {
         return (*v_)[n];
    }
    else
    {
         return (*delegate_).get(n+start_);
    }
}

template<class T>
inline void Foam::gpuList<T>::set(const label n, const T val)
{
    if(v_)
    {
        (*v_)[n] = val;
    }
    else
    {
        (*delegate_).set(n+start_,val);
    }
}



template<class T>
inline void Foam::gpuList<T>::clear()
{
    if(v_)
    {
         v_->clear();
    }
}

template<class T>
inline Foam::label Foam::gpuList<T>::size() const
{
    if(v_)
    {
        return v_->size();
    }
    else
    {
        return size_;
    }
}

template<class T>
inline void Foam::gpuList<T>::setSize(label size,const T val)
{
    if(v_)
    {
        v_->resize(size,val);
    }
    else
    {
         size_ = size;
    }
}

template<class T>
inline void Foam::gpuList<T>::setSize(label size)
{
    if(v_)
    {
        v_->resize(size);
    }
    else
    {
         size_ = size;
    }
}


template<class T>
inline bool Foam::gpuList<T>::empty() const
{
    return ! size();
}

