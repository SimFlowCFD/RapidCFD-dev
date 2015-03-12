#include "error.H"

#include "gpuList.H"
#include "contiguous.H"


template<class T>
Foam::gpuList<T>::~gpuList()
{
    if (this->v_)
    { 
        try
        {
             delete this->v_;
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
std::streamsize Foam::gpuList<T>::byteSize() const
{
    if (!contiguous<T>())
    {
        FatalErrorIn("gpuList<T>::byteSize()")
            << "Cannot return the binary size of a list of "
               "non-primitive elements"
            << abort(FatalError);
    }

    return this->size()*sizeof(T);
}

template<class T>
template<class Iterator>
void Foam::gpuList<T>::copyInto(Iterator it) const
{
    gpu_api::copy(this->begin(),this->end(),it);
}

template<class T>
void Foam::gpuList<T>::operator=(const T& t)
{
    gpu_api::fill(begin(),end(),t);
}

template<class T>
void Foam::gpuList<T>::operator=(const gpuList<T>& l)
{
    if(size() != l.size())
    {
        setSize(l.size());
    }

    gpu_api::copy(l.begin(),l.end(),begin());
}

template<class T>
void Foam::gpuList<T>::operator=(const UList<T>& l)
{
    setSize(l.size());
    gpu_api::copy(l.begin(),l.end(),begin());
}

template<class T>
void Foam::gpuList<T>::operator=(const List<T>& l)
{
    this->operator=(static_cast<const UList<T>&>(l));
}


template<class T>
void Foam::sort(gpuList<T>& a)
{
    gpu_api::sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::sort(gpuList<T>& a, const Cmp& cmp)
{
    gpu_api::sort(a.begin(), a.end(), cmp);
}


template<class T>
void Foam::stableSort(gpuList<T>& a)
{
    gpu_api::stable_sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::stableSort(gpuList<T>& a, const Cmp& cmp)
{
    gpu_api::stable_sort(a.begin(), a.end(), cmp);
}


#include "gpuListIO.C"

