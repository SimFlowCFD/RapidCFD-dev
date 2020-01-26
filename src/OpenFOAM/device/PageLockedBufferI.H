
template<class T>
inline void Foam::makePageLocked(UList<T>& buff, const label size)
{
    // this is a bit complicated since we do not have custom allocators
    T* mem = allocPageLocked<T>(size);
    UList<T> placeholder(mem, size);
    if(buff.size()) delete[] buff.data();
    transferData(placeholder, buff);
}

template<class T>
inline void Foam::clearPageLocked(UList<T>& buff)
{
    // UList does not free memory so we move it back
    UList<T> placeholder;
    transferData(buff, placeholder);
    freePageLocked<T>(placeholder.data());
}

template<class T>
inline Foam::PageLockedBuffer<T>::PageLockedBuffer():
    buffer_()
{}


template<class T>
inline Foam::PageLockedBuffer<T>::PageLockedBuffer(const label size):
    buffer_()
{
    makePageLocked<T>(buffer_, size);
}

template<class T>
inline Foam::PageLockedBuffer<T>::~PageLockedBuffer()
{
    if(buffer_.data()) clearPageLocked<T>(buffer_);
}

template<class T>
inline Foam::Field<T>& Foam::PageLockedBuffer<T>::buffer(const label size) const
{
    if (buffer_.size() != size)
    {
        if (buffer_.size()) 
            clearPageLocked<T>(buffer_);
        makePageLocked<T>(buffer_, size);
    }
    return buffer_;
}
