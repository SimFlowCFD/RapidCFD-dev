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

Class
    Foam::gpuField

Description
    Generic templated field type.

SourceFiles
    gpuFieldFunctions.H
    gpuFieldFunctionsM.H
    gpuFieldMapper.H
    gpuFieldM.H
    gpuField.C
    gpuFieldFunctions.C
    gpuFieldFunctionsM.C

\*---------------------------------------------------------------------------*/

#ifndef gpuField_H
#define gpuField_H

#include "tmp.H"
#include "direction.H"
#include "VectorSpace.H"
#include "gpuList.H"
#include "scalarList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class gpuField;

template<class Type>
class Field;

template<class Type>
Ostream& operator<<(Ostream&, const gpuField<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<gpuField<Type> >&);

class gpuFieldMapper;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class gpuField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gpuField
:
    public refCount,
    public gpuList<Type>
{

public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // Static data members

        static const char* const typeName;


    // Static Member Functions

        //- Return a null field
        inline static const gpuField<Type>& null()
        {
            return *reinterpret_cast< gpuField<Type>* >(0);
        }


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        gpuField();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit gpuField(const label);

        //- Construct given size and initial value
        explicit gpuField(const label, const Type&);

        //- Construct as copy of a gpuList\<Type\>
        explicit gpuField(const gpuList<Type>&);

        explicit gpuField(const Field<Type>&);

        //- Construct by transferring the List contents
        explicit gpuField(const Xfer<gpuList<Type> >&);

        //- Construct by 1 to 1 mapping from the given field
        gpuField
        (
            const gpuList<Type>& mapF,
            const labelgpuList& mapAddressing
        );

        //- Construct by 1 to 1 mapping from the given tmp field
        gpuField
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        gpuField
        (
            const gpuList<Type>& mapF,
            const labelgpuListList& mapAddressing,
            const scalargpuListList& weights
        );

        //- Construct by interpolative mapping from the given tmp field
        gpuField
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuListList& mapAddressing,
            const scalargpuListList& weights
        );

        //- Construct by mapping from the given field
        gpuField
        (
            const gpuList<Type>& mapF,
            const gpuFieldMapper& map
        );

        //- Construct by mapping from the given field
        gpuField
        (
            const gpuList<Type>& mapF,
            const gpuFieldMapper& map,
            const Type& defaultValue
        );

        //- Construct by mapping from the given field
        gpuField
        (
            const gpuList<Type>& mapF,
            const gpuFieldMapper& map,
            const gpuList<Type>& defaultValues
        );

        //- Construct by mapping from the given tmp field
        gpuField
        (
            const tmp<gpuField<Type> >& tmapF,
            const gpuFieldMapper& map
        );

        //- Construct by mapping from the given tmp field. Supplied uniform
        //  value for unmapped items
        gpuField
        (
            const tmp<gpuField<Type> >& tmapF,
            const gpuFieldMapper& map,
            const Type& defaultValue
        );

        //- Construct by mapping from the given tmp field. Supplied values
        //  for unmapped items
        gpuField
        (
            const tmp<gpuField<Type> >& tmapF,
            const gpuFieldMapper& map,
            const gpuList<Type>& defaultValues
        );

        //- Construct as copy
        gpuField(const gpuField<Type>&);

        //- Construct as copy or re-use as specified.
        explicit gpuField(gpuField<Type>&, bool reUse);

        //- Construct as sublist.
        explicit gpuField(const gpuList<Type>&, label size);

        //- Construct as sublist.
        explicit gpuField(const gpuList<Type>&, label size, label start);

        //- Construct by transferring the gpuField contents
        gpuField(const Xfer<gpuField<Type> >&);

        //- Construct as copy of tmp<gpuField>
#       ifdef ConstructFromTmp
        gpuField(const tmp<gpuField<Type> >&);
#       endif

        //- Construct from Istream
        gpuField(Istream&);

        //- Construct from a dictionary entry
        gpuField(const word& keyword, const dictionary&, const label size);

        //- Clone
        tmp<gpuField<Type> > clone() const;

        //- Return a pointer to a new calculatedFvPatchgpuFieldgpuField created on
        //  freestore without setting patchgpuField values
        template<class Type2>
        static tmp<gpuField<Type> > NewCalculatedType(const gpuField<Type2>& f)
        {
            return tmp<gpuField<Type> >(new gpuField<Type>(f.size()));
        }


    // Member Functions

        //- 1 to 1 map from the given field
        void map
        (
            const gpuList<Type>& mapF,
            const labelgpuList& mapAddressing
        );

        //- 1 to 1 map from the given tmp field
        void map
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuList& mapAddressing
        );

        //- Interpolative map from the given field
        void map
        (
            const gpuList<Type>& mapF,
            const labelgpuListList& mapAddressing,
            const scalargpuListList& weights
        );

        //- Interpolative map from the given tmp field
        void map
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuListList& mapAddressing,
            const scalargpuListList& weights
        );

        //- Map from the given field
        void map
        (
            const gpuList<Type>& mapF,
            const gpuFieldMapper& map
        );

        //- Map from the given tmp field
        void map
        (
            const tmp<gpuField<Type> >& tmapF,
            const gpuFieldMapper& map
        );

        //- Map from self
        void autoMap
        (
            const gpuFieldMapper& map
        );

        //- 1 to 1 reverse-map from the given field
        void rmap
        (
            const gpuList<Type>& mapF,
            const labelgpuList& mapAddressing
        );

        //- 1 to 1 reverse-map from the given tmp field
        void rmap
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuList& mapAddressing
        );

        //- Interpolative reverse map from the given field
        void rmap
        (
            const gpuList<Type>& mapF,
            const labelgpuList& mapAddressing,
            const gpuList<scalar>& weights
        );

        //- Interpolative reverse map from the given tmp field
        void rmap
        (
            const tmp<gpuField<Type> >& tmapF,
            const labelgpuList& mapAddressing,
            const gpuList<scalar>& weights
        );

        //- Negate this field
        void negate();

        //- Return a component field of the field
        tmp<gpuField<cmptType> > component(const direction) const;

        //- Replace a component field of the field
        void replace(const direction, const gpuList<cmptType>&);

        //- Replace a component field of the field
        void replace(const direction, const tmp<gpuField<cmptType> >&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        //- Return the field transpose (only defined for second rank tensors)
        tmp<gpuField<Type> > T() const;

        //- Write the field as a dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;

        tmp<Field<Type> > asField() const;


    // Member operators

        void operator=(const gpuField<Type>&);
        void operator=(const UList<Type>&);
	void operator=(const gpuList<Type>&);
        void operator=(const tmp<gpuField<Type> >&);
        void operator=(const Type&);

        template<class Form, class Cmpt, int nCmpt>
        void operator=(const VectorSpace<Form,Cmpt,nCmpt>&);

	void operator+=(const gpuList<Type>&);
        void operator+=(const tmp<gpuField<Type> >&);

	void operator-=(const gpuList<Type>&);
        void operator-=(const tmp<gpuField<Type> >&);

        void operator*=(const gpuList<scalar>&);
        void operator*=(const tmp<gpuField<scalar> >&);

        void operator/=(const gpuList<scalar>&);
        void operator/=(const tmp<gpuField<scalar> >&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);


    // IOstream operators

        friend Ostream& operator<< <Type>
        (Ostream&, const gpuField<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const tmp<gpuField<Type> >&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gpuFieldFunctions.H"

#ifdef NoRepository
#   include "gpuField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
