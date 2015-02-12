/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::UIndirectList

Description
    A List with indirect addressing.

    Like IndirectList but does not store addressing.

    Note the const_cast of the completeList. This is so we can use it both
    on const and non-const lists. Alternative would be to have a const_
    variant etc.

SourceFiles
    UIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef gpuIndirectList_H
#define gpuIndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T> class gpuIndirectList;
template<class T> Ostream& operator<<(Ostream&, const gpuIndirectList<T>&);

/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class gpuIndirectList
{
    // Private data

        gpuList<T>& completeList_;
        const labelgpuList& addressing_;


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline gpuIndirectList(const gpuList<T>&, const labelgpuList&);


    // Member Functions

        // Access

            //- Return the number of elements in the list
            inline label size() const;

            //- Return true if the list is empty (ie, size() is zero).
            inline bool empty() const;

            //- Return the complete list
            inline const gpuList<T>& completeList() const;

            //- Return the list addressing
            inline const gpuList<label>& addressing() const;


        // Member Operators
/*
            //- Return the addressed elements as a List
            inline List<T> operator()() const;

            //- Return non-const access to an element
            inline T& operator[](const label);

            //- Return const access to an element
            inline const T& operator[](const label) const;
*/
            //- Assignment from UList of addressed elements
            inline void operator=(const gpuList<T>&);

            //- Assignment from UIndirectList of addressed elements
            inline void operator=(const gpuIndirectList<T>&);

            //- Assignment of all entries to the given value
            inline void operator=(const T&);


    // STL type definitions

        //- Type of values the UList contains.
        typedef typename thrust::permutation_iterator<typename gpuList<T>::iterator,typename gpuList<label>::iterator> iterator;


        inline const iterator begin();
        inline const iterator begin() const;

        inline const iterator end();
        inline const iterator end() const;


    // Ostream operator

        //- Write UIndirectList to Ostream
        //  Binary output is currently still a bit of a problem
        friend Ostream& operator<< <T>
        (
            Ostream&,
            const gpuIndirectList<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gpuIndirectListI.H"

#ifdef NoRepository
#   include "gpuIndirectListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
