#pragma once

#include "gpuField.H"

namespace Foam
{

namespace cache
{


template<class Type>
inline gpuField<Type>& retrieve
(
    PtrList<gpuField<Type> >& list, 
    label level,
    label size
)
{
    if(level >= list.size())
        list.setSize(level+1);

    if(list.set(level))
    {
        gpuField<Type>& out = list[level];
        if(out.size() < size)
            out.setSize(size);
        return out; 
    }
    else
    {
        list.set(level,new gpuField<Type>(size));
        return list[level];
    }
}

template<class Type>
inline const gpuField<Type>& retrieveConst
(
    PtrList<gpuField<Type> >& list, 
    label level,
    label size
)
{
    return const_cast<const gpuField<Type>&>(retrieve(list,level,size));
}


}
}
