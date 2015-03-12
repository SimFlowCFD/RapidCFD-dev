/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gpuFieldMapper.H"
#include "gpuFieldM.H"
#include "dictionary.H"
#include "contiguous.H"
#include "gpuField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class Type>
const char* const Foam::gpuField<Type>::typeName("gpuField");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::gpuField<Type>::gpuField()
:
    gpuList<Type>()
{}


template<class Type>
Foam::gpuField<Type>::gpuField(const label size)
:
    gpuList<Type>(size)
{}


template<class Type>
Foam::gpuField<Type>::gpuField(const label size, const Type& t)
:
    gpuList<Type>(size, t)
{}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const gpuList<Type>& mapF,
    const labelgpuList& mapAddressing
)
:
    gpuList<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuList& mapAddressing
)
:
    gpuList<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const gpuList<Type>& mapF,
    const labelgpuListList& mapAddressing,
    const scalargpuListList& mapWeights
)
:
    gpuList<Type>(mapAddressing.size())
{
    map(mapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuListList& mapAddressing,
    const scalargpuListList& mapWeights
)
:
    gpuList<Type>(mapAddressing.size())
{
    map(tmapF, mapAddressing, mapWeights);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const gpuList<Type>& mapF,
    const gpuFieldMapper& mapper
)
:
    gpuList<Type>(mapper.size())
{
    map(mapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const gpuList<Type>& mapF,
    const gpuFieldMapper& mapper,
    const Type& defaultValue
)
:
    gpuList<Type>(mapper.size(), defaultValue)
{
    map(mapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const gpuList<Type>& mapF,
    const gpuFieldMapper& mapper,
    const gpuList<Type>& defaultValues
)
:
    gpuList<Type>(defaultValues)
{
    map(mapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const tmp<gpuField<Type> >& tmapF,
    const gpuFieldMapper& mapper
)
:
    gpuList<Type>(mapper.size())
{
    map(tmapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const tmp<gpuField<Type> >& tmapF,
    const gpuFieldMapper& mapper,
    const Type& defaultValue
)
:
    gpuList<Type>(mapper.size(), defaultValue)
{
    map(tmapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const tmp<gpuField<Type> >& tmapF,
    const gpuFieldMapper& mapper,
    const gpuList<Type>& defaultValues
)
:
    gpuList<Type>(defaultValues)
{
    map(tmapF, mapper);
}


template<class Type>
Foam::gpuField<Type>::gpuField(const gpuField<Type>& f)
:
    refCount(),
    gpuList<Type>(f)
{
}

template<class Type>
Foam::gpuField<Type>::gpuField(const Field<Type>& f)
:
    gpuList<Type>(f.size())
{
    this->operator=(f);
}


template<class Type>
Foam::gpuField<Type>::gpuField(gpuField<Type>& f, bool reUse)
:
    gpuList<Type>(static_cast<gpuList<Type>&>(f), reUse)
{}

template<class Type>
Foam::gpuField<Type>::gpuField(const gpuField<Type>& f, label size)
:
    gpuList<Type>(static_cast<const gpuList<Type>&>(f), size)
{}

template<class Type>
Foam::gpuField<Type>::gpuField(const gpuField<Type>& f, label size, label start)
:
    gpuList<Type>(static_cast<const gpuList<Type>&>(f), size, start)
{}


template<class Type>
Foam::gpuField<Type>::gpuField(const gpuList<Type>& l, label size)
:
    gpuList<Type>(l, size)
{}

template<class Type>
Foam::gpuField<Type>::gpuField(const gpuList<Type>& l, label size, label start)
:
    gpuList<Type>(l, size, start)
{}


template<class Type>
Foam::gpuField<Type>::gpuField(const Xfer<gpuList<Type> >& f)
:
    gpuList<Type>(f)
{}


template<class Type>
Foam::gpuField<Type>::gpuField(const Xfer<gpuField<Type> >& f)
:
    gpuList<Type>()
{
    gpuList<Type>::transfer(f());	
}


template<class Type>
Foam::gpuField<Type>::gpuField(const gpuList<Type>& list)
:
    gpuList<Type>(list)
{}


// Construct as copy of tmp<gpuField>
#ifdef ConstructFromTmp
template<class Type>
Foam::gpuField<Type>::gpuField(const tmp<gpuField<Type> >& tf)
:
    gpuList<Type>(const_cast<gpuField<Type>&>(tf()), tf.isTmp())
{
    const_cast<gpuField<Type>&>(tf()).resetRefCount();
}
#endif


template<class Type>
Foam::gpuField<Type>::gpuField(Istream& is)
:
    gpuList<Type>(is)
{}


template<class Type>
Foam::gpuField<Type>::gpuField
(
    const word& keyword,
    const dictionary& dict,
    const label s
)
{
    Field<Type> tf(keyword,dict,s);
    this->operator=(tf);
}


template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::gpuField<Type>::clone() const
{
    return tmp<gpuField<Type> >(new gpuField<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::gpuField<Type>::map
(
    const gpuList<Type>& mapF,
    const labelgpuList& mapAddressing
)
{
    gpuField<Type>& f = *this;

    if (f.size() != mapAddressing.size())
    {
        f.setSize(mapAddressing.size());
    }

    
    if (mapF.size() > 0)
    {
        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                mapF.begin(), 
                mapAddressing.begin()
            ),
            thrust::make_permutation_iterator
            (
                mapF.begin(), 
                mapAddressing.end()
            ),
            f.begin()
        );
    }
}


template<class Type>
void Foam::gpuField<Type>::map
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuList& mapAddressing
)
{
    map(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::gpuField<Type>::map
(
    const gpuList<Type>& mapF,
    const labelgpuListList& mapAddressing,
    const scalargpuListList& mapWeights
)
{
    gpuField<Type>& f = *this;

    if (f.size() != mapAddressing.size())
    {
        f.setSize(mapAddressing.size());
    }

    if (mapWeights.size() != mapAddressing.size())
    {
        FatalErrorIn
        (
            "void gpuField<Type>::map\n"
            "(\n"
            "    const gpuList<Type>& mapF,\n"
            "    const labelgpuListList& mapAddressing,\n"
            "    const scalargpuListList& mapWeights\n"
            ")"
        ) << "Weights and addressing map have different sizes.  Weights size: "
            << mapWeights.size() << " map size: " << mapAddressing.size()
            << abort(FatalError);
    }

    //TODO transform into gpu code
    notImplemented
    (
        "Type Foam::gpuField<Type>::map"
        "("
             "const gpuList<Type>& mapF,"
             "const labelgpuListList& mapAddressing,"
             "const scalargpuListList& mapWeights"
        ")"
    );
/*
    forAll(f, i)
    {
        const labelgpuList&  localAddrs   = mapAddressing[i];
        const scalargpuList& localWeights = mapWeights[i];

        f[i] = pTraits<Type>::zero;

        forAll(localAddrs, j)
        {
            f[i] += localWeights[j]*mapF[localAddrs[j]];
        }
    }
*/
}


template<class Type>
void Foam::gpuField<Type>::map
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuListList& mapAddressing,
    const scalargpuListList& mapWeights
)
{
    map(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::gpuField<Type>::map
(
    const gpuList<Type>& mapF,
    const gpuFieldMapper& mapper
)
{
    if
    (
        mapper.direct()
     && &mapper.directAddressing()
     && mapper.directAddressing().size()
    )
    {
        map(mapF, mapper.directAddressing());
    }
    else if (!mapper.direct() && mapper.addressing().size())
    {
        map(mapF, mapper.addressing(), mapper.weights());
    }
}


template<class Type>
void Foam::gpuField<Type>::map
(
    const tmp<gpuField<Type> >& tmapF,
    const gpuFieldMapper& mapper
)
{
    map(tmapF(), mapper);
    tmapF.clear();
}


template<class Type>
void Foam::gpuField<Type>::autoMap
(
    const gpuFieldMapper& mapper
)
{
    if
    (
        (
            mapper.direct()
         && &mapper.directAddressing()
         && mapper.directAddressing().size()
        )
     || (!mapper.direct() && mapper.addressing().size())
    )
    {
        gpuField<Type> fCpy(*this);
        map(fCpy, mapper);
    }
    else
    {
        this->setSize(mapper.size());
    }
}


template<class Type>
void Foam::gpuField<Type>::rmap
(
    const gpuList<Type>& mapF,
    const labelgpuList& mapAddressing
)
{
    gpuField<Type>& f = *this;
    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            mapF.begin(), 
            mapAddressing.begin()
        ),
        thrust::make_permutation_iterator
        (
            mapF.begin(), 
            mapAddressing.end()
        ),
        f.begin()
    );
}


template<class Type>
void Foam::gpuField<Type>::rmap
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuList& mapAddressing
)
{
    rmap(tmapF(), mapAddressing);
    tmapF.clear();
}


template<class Type>
void Foam::gpuField<Type>::rmap
(
    const gpuList<Type>& mapF,
    const labelgpuList& mapAddressing,
    const gpuList<scalar>& mapWeights
)
{
    gpuField<Type>& f = *this;
    f = pTraits<Type>::zero;
    gpuList<Type> tmp(f.size());
    tmp = pTraits<Type>::zero;

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            mapF.begin(), 
            mapAddressing.begin()
        ),
        thrust::make_permutation_iterator
        (
            mapF.begin(), 
            mapAddressing.end()
        ),
        tmp.begin()
    );

    thrust::transform
    (
        tmp.begin(),
        tmp.end(),
        mapWeights.begin(),
        f.begin(),
        multiplyOperatorFunctor<Type,scalar,Type>()
    );

}


template<class Type>
void Foam::gpuField<Type>::rmap
(
    const tmp<gpuField<Type> >& tmapF,
    const labelgpuList& mapAddressing,
    const gpuList<scalar>& mapWeights
)
{
    rmap(tmapF(), mapAddressing, mapWeights);
    tmapF.clear();
}


template<class Type>
void Foam::gpuField<Type>::negate()
{
    thrust::transform(this->begin(),this->end(),this->begin(),thrust::negate<Type>());
}


template<class Type>
Foam::tmp<Foam::gpuField<typename Foam::gpuField<Type>::cmptType> >
Foam::gpuField<Type>::component
(
    const direction d
) const
{
    tmp<gpuField<cmptType> > Component(new gpuField<cmptType>(this->size()));
    ::Foam::component(Component(), *this, d);
    return Component;
}


template<class Type>
void Foam::gpuField<Type>::replace
(
    const direction d,
    const gpuList<cmptType>& sf
)
{
    thrust::transform(this->begin(),this->end(),sf.begin(),this->begin(),
            replaceComponentFunctor<Type,cmptType>(d));
}


template<class Type>
void Foam::gpuField<Type>::replace
(
    const direction d,
    const tmp<gpuField<cmptType> >& tsf
)
{
    replace(d, tsf());
    tsf.clear();
}


template<class Type>
void Foam::gpuField<Type>::replace
(
    const direction d,
    const cmptType& c
)
{
    thrust::transform(this->begin(),this->end(),this->begin(),
              replaceComponentWithSourceFunctor<Type,cmptType>(d,c));
}


template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::gpuField<Type>::T() const
{
    tmp<gpuField<Type> > transpose(new gpuField<Type>(this->size()));
    ::Foam::T(transpose(), *this);
    return transpose;
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::gpuField<Type>::asField() const
{
    Field<Type>* f = new Field<Type>(this->size());

    thrust::copy(this->begin(),this->end(),f->begin());

    return tmp<Field<Type> >(f);
}


template<class Type>
void Foam::gpuField<Type>::writeEntry(const word& keyword, Ostream& os) const
{
    Field<Type> f(this->asField());
    
    f.writeEntry(keyword,os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::gpuField<Type>::operator=(const gpuField<Type>& rhs)
{
    if (this == &rhs)
    {
        FatalErrorIn("gpuField<Type>::operator=(const gpuField<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    gpuList<Type>::operator=(rhs);
}

template<class Type>
void Foam::gpuField<Type>::operator=(const gpuList<Type>& rhs)
{
    gpuList<Type>::operator=(rhs);
}

template<class Type>
void Foam::gpuField<Type>::operator=(const UList<Type>& rhs)
{
    gpuList<Type>::operator=(rhs);
}

template<class Type>
void Foam::gpuField<Type>::operator=(const tmp<gpuField>& rhs)
{
    if (this == &(rhs()))
    {
        FatalErrorIn("gpuField<Type>::operator=(const tmp<gpuField>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    // This is dodgy stuff, don't try it at home.
    gpuField* fieldPtr = rhs.ptr();
    gpuList<Type>::transfer(*fieldPtr);
    delete fieldPtr;
}


template<class Type>
void Foam::gpuField<Type>::operator=(const Type& t)
{
    gpuList<Type>::operator=(t);
}


template<class Type>
template<class Form, class Cmpt, int nCmpt>
void Foam::gpuField<Type>::operator=(const VectorSpace<Form,Cmpt,nCmpt>& vs)
{
    thrust::fill(this->begin(),this->end(),vs);
}


#define COMPUTED_ASSIGNMENT(TYPE, op, opFunc)                                 \
                                                                              \
template<class Type>                                                          \
void Foam::gpuField<Type>::operator op(const gpuList<TYPE>& f)                \
{                                                                             \
    thrust::transform(this->begin(), this->end(), f.begin(), this->begin(),   \
                   opFunc##OperatorFunctor<Type,TYPE,Type>());                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::gpuField<Type>::operator op(const tmp<gpuField<TYPE> >& tf)        \
{                                                                             \
    operator op(tf());                                                        \
    tf.clear();                                                               \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void Foam::gpuField<Type>::operator op(const TYPE& t)                         \
{                                                                             \
    thrust::transform(this->begin(), this->end(), this->begin(),              \
                   opFunc##OperatorFSFunctor<Type,TYPE,Type>(t));             \
}

COMPUTED_ASSIGNMENT(Type, +=, add)
COMPUTED_ASSIGNMENT(Type, -=, subtract)
COMPUTED_ASSIGNMENT(scalar, *=, multiply)
COMPUTED_ASSIGNMENT(scalar, /=, divide)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const gpuField<Type>& f)
{
    os  << static_cast<const gpuList<Type>&>(f);
    return os;
}


template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const tmp<gpuField<Type> >& tf)
{
    Field<Type> f(tf().asField());
    os  << f;
    tf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gpuFieldFunctions.C"

// ************************************************************************* //
