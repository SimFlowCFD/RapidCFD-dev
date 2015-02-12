/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::edge

Description
    An edge is a list of two point labels. The functionality it provides
    supports the discretisation on a 2-D flat mesh.

SourceFiles
    edgeI.H

\*---------------------------------------------------------------------------*/

#ifndef edge_H
#define edge_H

#include "FixedList.H"
#include "pointField.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class edge;
__HOST____DEVICE__
inline bool operator==(const edge& a, const edge& b);
__HOST____DEVICE__
inline bool operator!=(const edge& a, const edge& b);

/*---------------------------------------------------------------------------*\
                           Class edge Declaration
\*---------------------------------------------------------------------------*/

class edge
:
    public FixedList<label, 2>
{

public:

    // Static data members

        static const char* const typeName;
        static const edge zero;


    // Constructors

        //- Null constructor for lists
        __HOST____DEVICE__
        inline edge();

        //- Construct from components
        __HOST____DEVICE__
        inline edge(const label a, const label b);

        //- Construct from FixedList
        __HOST____DEVICE__
        inline edge(const FixedList<label, 2>&);

        //- Construct from Istream
        inline edge(Istream&);


    // Member Functions

        //- Return start vertex label
        __HOST____DEVICE__
        inline label start() const;

        //- Return start vertex label
        __HOST____DEVICE__
        inline label& start();

        //- Return end vertex label
        __HOST____DEVICE__
        inline label end() const;

        //- Return end vertex label
        __HOST____DEVICE__
        inline label& end();

        //- Given one vertex, return the other
        __HOST____DEVICE__
        inline label otherVertex(const label a) const;

        //- Return common vertex
        __HOST____DEVICE__
        inline label commonVertex(const edge& a) const;

        //- Flip the edge in-place.
        __HOST____DEVICE__
        inline void flip();

        //- Return reverse edge
        __HOST____DEVICE__
        inline edge reverseEdge() const;

        //- Return centre (centroid)
        __HOST____DEVICE__
        inline point centre(const point*) const;
        inline point centre(const pointField&) const;

        //- Return the vector (end - start)
        __HOST____DEVICE__
        inline vector vec(const point*) const;
        inline vector vec(const pointField&) const;

        //- Return scalar magnitude
        __HOST____DEVICE__
        inline scalar mag(const point*) const;
        inline scalar mag(const pointField&) const;

        //- Return edge line
        __HOST____DEVICE__
        inline linePointRef line(const point*) const;
        inline linePointRef line(const pointField&) const;

        //- compare edges
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        __HOST____DEVICE__
        static inline int compare(const edge&, const edge&);


    // Friend Operators

        __HOST____DEVICE__
        friend bool operator==(const edge& a, const edge& b);
        __HOST____DEVICE__
        friend bool operator!=(const edge& a, const edge& b);
};


//- Hash specialization for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e, unsigned seed) const
{
    unsigned val = seed;

    if (e.start() < e.end())
    {
        val = Hash<label>()(e.start(), val);
        val = Hash<label>()(e.end(), val);
    }
    else
    {
        val = Hash<label>()(e.end(), val);
        val = Hash<label>()(e.start(), val);
    }

    return val;
}


//- Hash specialization for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e) const
{
    return Hash<edge>()(e, 0);
}


template<>
inline bool contiguous<edge>()  {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
