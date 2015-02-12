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
    Foam::quaternion

Description
    Quaternion class used to perform rotations in 3D space.

SourceFiles
    quaternionI.H
    quaternion.C

\*---------------------------------------------------------------------------*/

#ifndef quaternion_H
#define quaternion_H

#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "word.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class quaternion;
Istream& operator>>(Istream& is, quaternion&);
Ostream& operator<<(Ostream& os, const quaternion& C);


/*---------------------------------------------------------------------------*\
                           Class quaternion Declaration
\*---------------------------------------------------------------------------*/

class quaternion
{
    // private data

        //- Scalar part of the quaternion ( = cos(theta/2) for rotation)
        scalar w_;

        //- Vector part of the quaternion ( = axis of rotation)
        vector v_;


        //- Multiply vector v by quaternion as if v is a pure quaternion
        __HOST____DEVICE__
        inline quaternion mulq0v(const vector& v) const;


public:

    //- Component type
    typedef scalar cmptType;

    // Member constants

        enum
        {
            rank = 1 // Rank of quaternion is 1
        };


    // Static data members

        static const char* const typeName;

        static const quaternion zero;
        static const quaternion I;


    // Constructors

        //- Construct null
        __HOST____DEVICE__
        inline quaternion();

        //- Construct given scalar and vector parts
        __HOST____DEVICE__
        inline quaternion(const scalar w, const vector& v);

        //- Construct a rotation quaternion given the direction d
        //  and angle theta
        __HOST____DEVICE__
        inline quaternion(const vector& d, const scalar theta);

        //- Construct a rotation quaternion given the direction d
        //  and cosine angle cosTheta and a if d is normalized
        __HOST____DEVICE__
        inline quaternion
        (
            const vector& d,
            const scalar cosTheta,
            const bool normalized
        );

        //- Construct given scalar part, the vector part = vector::zero
        __HOST____DEVICE__
        inline explicit quaternion(const scalar w);

        //- Construct a pure quaternion given the vector part, scalar part = 0
        __HOST____DEVICE__
        inline explicit quaternion(const vector& v);

        //- Construct a quaternion given the three Euler angles
        __HOST____DEVICE__
        inline quaternion
        (
            const scalar angleX,
            const scalar angleY,
            const scalar angleZ
        );

        //- Construct a quaternion from a rotation tensor
        __HOST____DEVICE__
        inline explicit quaternion(const tensor& rotationTensor);

        //- Construct from Istream
        quaternion(Istream&);


    // Member functions

           // Access

               //- Scalar part of the quaternion ( = cos(theta/2) for rotation)
               __HOST____DEVICE__
               inline scalar w() const;

               //- Vector part of the quaternion ( = axis of rotation)
               __HOST____DEVICE__
               inline const vector& v() const;

               //- The rotation tensor corresponding the quaternion
               __HOST____DEVICE__
               inline tensor R() const;

               //- Return a vector of euler angles (rotations in radians about
               //  the x, y and z axes.
               __HOST____DEVICE__
               inline vector eulerAngles(const quaternion& q) const;

               __HOST____DEVICE__
               inline quaternion normalized() const;


           // Edit

               //- Scalar part of the quaternion ( = cos(theta/2) for rotation)
               __HOST____DEVICE__
               inline scalar& w();

               //- Vector part of the quaternion ( = axis of rotation)
               __HOST____DEVICE__
               inline vector& v();

               __HOST____DEVICE__
               inline void normalize();


           // Transform

               //- Rotate the given vector
               __HOST____DEVICE__
               inline vector transform(const vector& v) const;

               //- Rotate the given vector anti-clockwise
               __HOST____DEVICE__
               inline vector invTransform(const vector& v) const;

               //- Rotate the given quaternion (and normalize)
               __HOST____DEVICE__
               inline quaternion transform(const quaternion& q) const;

               //- Rotate the given quaternion anti-clockwise (and normalize)
               __HOST____DEVICE__
               inline quaternion invTransform(const quaternion& q) const;


    // Member operators

        __HOST____DEVICE__
        inline void operator=(const quaternion&);
        __HOST____DEVICE__
        inline void operator+=(const quaternion&);
        __HOST____DEVICE__
        inline void operator-=(const quaternion&);
        __HOST____DEVICE__
        inline void operator*=(const quaternion&);
        __HOST____DEVICE__
        inline void operator/=(const quaternion&);

        __HOST____DEVICE__
        inline void operator=(const scalar);

        __HOST____DEVICE__
        inline void operator=(const vector&);

        __HOST____DEVICE__
        inline void operator*=(const scalar);
        __HOST____DEVICE__
        inline void operator/=(const scalar);


    // IOstream operators

        friend Istream& operator>>(Istream& is, quaternion&);
        friend Ostream& operator<<(Ostream& os, const quaternion& C);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline scalar magSqr(const quaternion& q);
__HOST____DEVICE__
inline scalar mag(const quaternion& q);

//- Return the conjugate of the given quaternion
__HOST____DEVICE__
inline quaternion conjugate(const quaternion& q);

//- Return the normalized (unit) quaternion of the given quaternion
__HOST____DEVICE__
inline quaternion normalize(const quaternion& q);

//- Return the inverse of the given quaternion
__HOST____DEVICE__
inline quaternion inv(const quaternion& q);

//- Return a string representation of a quaternion
word name(const quaternion&);

//- Spherical linear interpolation of quaternions
__HOST____DEVICE__
inline quaternion slerp
(
    const quaternion& qa,
    const quaternion& qb,
    const scalar t
);

//- Exponent of a quaternion
__HOST____DEVICE__
inline quaternion exp(const quaternion& q);

//- Power of a quaternion
__HOST____DEVICE__
inline quaternion pow(const quaternion& q, const label power);

//- Power of a quaternion
__HOST____DEVICE__
inline quaternion pow(const quaternion& q, const scalar power);

//- Data associated with quaternion type are contiguous
template<>
__HOST____DEVICE__
inline bool contiguous<quaternion>() {return true;}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline bool operator==(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline bool operator!=(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline quaternion operator+(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline quaternion operator-(const quaternion& q);
__HOST____DEVICE__
inline quaternion operator-(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline scalar operator&(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline quaternion operator*(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline quaternion operator/(const quaternion& q1, const quaternion& q2);
__HOST____DEVICE__
inline quaternion operator*(const scalar s, const quaternion& q);
__HOST____DEVICE__
inline quaternion operator*(const quaternion& q, const scalar s);
__HOST____DEVICE__
inline quaternion operator/(const quaternion& q, const scalar s);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quaternionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
