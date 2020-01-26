#pragma once

template<class T>
inline T* Foam::gpuList<T>::data()
{
    return v_;
}


template<class T>
inline const T* Foam::gpuList<T>::data() const
{
    return v_;
}


template<class T>
inline typename Foam::gpuList<T>::iterator Foam::gpuList<T>::begin()
{
    return iterator(v_);
}


template<class T>
inline typename Foam::gpuList<T>::iterator Foam::gpuList<T>::end()
{
    return iterator(v_+size_);
}


template<class T>
inline typename Foam::gpuList<T>::const_iterator Foam::gpuList<T>::begin() const
{
    return const_iterator(v_);
}


template<class T>
inline typename Foam::gpuList<T>::const_iterator Foam::gpuList<T>::end() const
{
    return const_iterator(v_+size_);
}


template<class T>
inline std::streamsize Foam::gpuList<T>::byteSize() const
{
    return size()*sizeof(T);
}


template<class T>
inline Foam::label Foam::gpuList<T>::size() const
{
    return size_;
}


template<class T>
inline bool Foam::gpuList<T>::empty() const
{
    return ! size_;
}
