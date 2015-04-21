#ifndef gpuList_H
#define gpuList_H

#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "Xfer.H"
#include "gpuConfig.H"

namespace Foam
{

template<class T> class List;
template<class T> class UList;
template<class T> class gpuList;

template<class T> Ostream& operator<<(Ostream&, const gpuList<T>&);
template<class T> Istream& operator>>(Istream&, gpuList<T>&);

typedef gpuList<label> labelgpuList;

template<class T>
class gpuList
{
        label size_;
        label start_;

        gpuList<T>* delegate_;
        gpu_api::device_vector<T>* v_;

public:

        inline static const gpuList<T>& null();


        inline gpuList();
        explicit inline gpuList(label size);
        explicit inline gpuList(label size, const T&);
        explicit inline gpuList(const Xfer<gpuList<T> >&);
        explicit inline gpuList(gpuList<T>& a, bool reUse);
        explicit inline gpuList(const gpuList<T>& list);
        explicit inline gpuList(const UList<T>& list);

        template<class InputIterator>
        inline gpuList(InputIterator first, InputIterator last);

        //- Construct sublist from gpuList and sub-list size
        explicit inline gpuList(const gpuList<T>& list,const label subSize);

        //- Construct sublist from gpuList, sub-list size and start index
        explicit inline gpuList
        (
            const gpuList<T>& list,
            const label subSize,
            const label startIndex
        );

        gpuList(Istream&);

        ~gpuList();


        void writeEntry(Ostream&) const;
        void writeEntry(const word& keyword, Ostream&) const;


        inline Xfer<gpuList<T> > xfer();
        void transfer(gpuList<T>&);

        void setDelegate(gpuList<T>&);
        void setDelegate(gpuList<T>&,label);
        void setDelegate(gpuList<T>&,label,label);

        //Requires iterator to have sufficient capacity
        template<class Iterator>
        void copyInto(Iterator it) const;

        void operator=(const T&);
        void operator=(const gpuList<T>&);
        void operator=(const UList<T>&);
        void operator=(const List<T>&);

        std::streamsize byteSize() const;
        inline T* data();
        inline const T* data() const;

        typedef typename gpu_api::device_vector<T>::iterator        iterator;
        typedef typename gpu_api::device_vector<T>::const_iterator        const_iterator;
        typedef typename gpu_api::device_vector<T>::reverse_iterator        reverse_iterator;
        typedef typename gpu_api::device_vector<T>::const_reverse_iterator        const_reverse_iterator;

        inline const iterator begin();
        inline const iterator end();

        inline const const_iterator begin() const;
        inline const const_iterator end() const;

        inline const reverse_iterator rbegin();
        inline const reverse_iterator rend();

        inline const const_reverse_iterator rbegin() const;
        inline const const_reverse_iterator rend() const;

        inline T first() const;

        inline T get(const label n) const;
        inline void set(const label n, const T val);

        inline void clear();
        inline label size() const;
        inline void setSize(label size);
        inline void setSize(label size, const T val);
        inline bool empty() const;


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

template<class T>
void sort(gpuList<T>&);

template<class T, class Cmp>
void sort(gpuList<T>&, const Cmp&);

template<class T>
void stableSort(gpuList<T>&);

template<class T, class Cmp>
void stableSort(gpuList<T>&, const Cmp&);

} 

#   include "gpuListI.H"

#ifdef NoRepository
#   include "gpuList.C"
#endif

#endif

// ************************************************************************* //
