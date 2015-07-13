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
    Foam::pointConstraint

Description
    Accumulates point constraints through successive applications of the
    applyConstraint function.

    After all the constraints have been entered the resulting
    transformation tensor is returned by the constraintTransformation
    function.

SourceFiles
    pointConstraintI.H

\*---------------------------------------------------------------------------*/

#ifndef pointConstraint_H
#define pointConstraint_H

#include "label.H"
#include "vector.H"
#include "transform.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointConstraint;

Istream& operator>>
(
    Istream& is,
    pointConstraint& pc
);

Ostream& operator<<
(
    Ostream& os,
    const pointConstraint& pc
);

/*---------------------------------------------------------------------------*\
                           Class pointConstraint Declaration
\*---------------------------------------------------------------------------*/

class pointConstraint
{
    label point_;
    vector value_;

public:

    // Constructors

        //- Construct null
        __host__ __device__
        inline pointConstraint();

        //- Construct from components
        __host__ __device__
        inline pointConstraint(const label&, const vector&);

        //- Construct from components
        inline pointConstraint(const Tuple2<label, vector>&);

        //- Construct from Istream
        inline pointConstraint(Istream&);


    // Member Functions

        //- Apply and accumulate the effect of the given constraint direction
        __host__ __device__
        inline void applyConstraint(const vector& cd);

        //- Combine constraints
        __host__ __device__
        inline void combine(const pointConstraint&);

        //- Return the accumulated constraint transformation tensor
        __host__ __device__
        inline tensor constraintTransformation() const;

        //- Return the accumulated unconstrained directions. Directions
        //  coded as first n rows of tensor.
        __host__ __device__
        inline void unconstrainedDirections(label& n, tensor& vecs) const;

        __host__ __device__
        inline label& first()
        {
            return point_;
        }

        __host__ __device__
        inline const label& first() const
        {
            return point_;
        }

        __host__ __device__
        inline vector& second()
        {
            return value_;
        }

        __host__ __device__
        inline const vector& second() const
        {
            return value_;
        }

    // IOstream operators

        //- Read Tuple2 from Istream, discarding contents of existing Tuple2.
        friend Istream& operator>>
        (
            Istream& is,
            pointConstraint& pc
        );

        // Write Tuple2 to Ostream.
        friend Ostream& operator<<
        (
            Ostream& os,
            const pointConstraint& pc
        );
};

inline bool operator==
(
    const pointConstraint& a,
    const pointConstraint& b
)
{
    return (a.first() == b.first() && a.second() == b.second());
}

inline bool operator!=
(
    const pointConstraint& a,
    const pointConstraint& b
)
{
    return !(a == b);
}

//- reduce operator
class combineConstraintsEqOp
{
public:
    __host__ __device__
    inline void operator()(pointConstraint&, const pointConstraint&) const;
};

//- transformation function
__host__ __device__
inline pointConstraint transform(const tensor& tt, const pointConstraint& v);

//- contiguous
template<class T> bool contiguous();
template<>
inline bool contiguous<pointConstraint>()                {return true;}


inline Istream& operator>>(Istream& is, pointConstraint& t2)
{
    is.readBegin("Tuple2");
    is >> t2.point_ >> t2.value_;
    is.readEnd("Tuple2");

    // Check state of Istream
    is.check("operator>>(Istream&, pointConstraint&)");

    return is;
}


inline Ostream& operator<<(Ostream& os, const pointConstraint& t2)
{
    os  << token::BEGIN_LIST
        << t2.point_ << token::SPACE << t2.value_
        << token::END_LIST;

    return os;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointConstraintI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
