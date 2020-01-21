#pragma once

#include "error.H"
#include "Xfer.H"
#include <thrust/device_ptr.h>

namespace Foam
{

template<class T> class List;
template<class T> class UList;
template<class T> class gpuList;

template<class T>
using device_ptr = thrust::device_ptr<T>;

template<class T> Ostream& operator<<(Ostream&, const gpuList<T>&);
template<class T> Istream& operator>>(Istream&, gpuList<T>&);

typedef gpuList<label> labelgpuList;

template<class T>
class gpuList
{
        T* v_;
        label size_;
        //- Am I the owner of this memory
        bool owner_;

public:

        static const gpuList<T>& null();


        gpuList();
        explicit gpuList(label size);
        explicit gpuList(label size, const T&);
        explicit gpuList(const Xfer<gpuList<T> >&);
        explicit gpuList(gpuList<T>& a, bool reUse);
        explicit gpuList(const gpuList<T>& list);
        explicit gpuList(const UList<T>& list);

        //- Construct sublist from gpuList and sub-list size
        explicit gpuList(const gpuList<T>& list,const label subSize);

        //- Construct sublist from gpuList, sub-list size and start index
        explicit gpuList
        (
            const gpuList<T>& list,
            const label subSize,
            const label startIndex
        );

        gpuList(Istream&);

        ~gpuList();


        void writeEntry(Ostream&) const;
        void writeEntry(const word& keyword, Ostream&) const;


        Xfer<gpuList<T>> xfer();
        void transfer(gpuList<T>&);

        void setDelegate(gpuList<T>&);
        void setDelegate(gpuList<T>&,label);
        void setDelegate(gpuList<T>&,label,label);

        void operator=(const T&);
        void operator=(const gpuList<T>&);
        void operator=(const UList<T>&);
        void operator=(const List<T>&);

        inline label size() const;
        inline bool empty() const;
        inline std::streamsize byteSize() const;
        inline T* data();
        inline const T* data() const;

        typedef device_ptr<T>        iterator;
        typedef device_ptr<const T>  const_iterator;

        inline iterator begin();
        inline iterator end();

        inline const_iterator begin() const;
        inline const_iterator end() const;

        T first() const;

        T get(const label n) const;
        void set(const label n, const T val);

        void clear();
        void setSize(label size);
        void setSize(label size, const T val);


        friend Ostream& operator<< <T>
        (
            Ostream&,
            const gpuList<T>&
        );

        friend Istream& operator>> <T>
        (
            Istream&,
            gpuList<T>&
        );
};

extern template class gpuList<bool>;
extern template class gpuList<char>;
extern template class gpuList<label>;
extern template class gpuList<float>;
extern template class gpuList<double>;

} 

#include "gpuListI.H"
