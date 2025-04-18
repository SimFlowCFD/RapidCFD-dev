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
    Foam::PointHit

Description
    This class describes the interaction of a face and a point. It
    carries the info of a successful hit and (if successful), returns
    the interaction point.

\*---------------------------------------------------------------------------*/

#ifndef PointHit_H
#define PointHit_H

#include "bool.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class Ostream;


// Forward declaration of friend functions and operators

template<class Point> class PointHit;

template<class Point>
inline Ostream& operator<<(Ostream&, const PointHit<Point>&);


/*---------------------------------------------------------------------------*\
                           Class PointHit Declaration
\*---------------------------------------------------------------------------*/

template<class Point>
class PointHit
{
    // Private data

        //- Hit success
        bool hit_;

        //- Point of hit; for miss holds best estimate outside the object
        Point hitPoint_;

        //- Distance to hit point
        scalar distance_;

        //- Eligible miss
        bool eligibleMiss_;


public:

    // Constructors

        //- Construct null
        __HOST____DEVICE__
        PointHit()
        :
            hit_(false),
            hitPoint_(vector(0.0,0.0,0.0)),
            distance_(GREAT),
            eligibleMiss_(false)
        {}

        //- Construct from components
        __HOST____DEVICE__
        PointHit
        (
            const bool hit,
            const Point& p,
            const scalar dist,
            const bool eligibleMiss
        )
        :
            hit_(hit),
            hitPoint_(p),
            distance_(dist),
            eligibleMiss_(eligibleMiss)
        {}

        //- Construct from point. Hit and distance set later
        __HOST____DEVICE__
        PointHit(const Point& p)
        :
            hit_(false),
            hitPoint_(p),
            distance_(GREAT),
            eligibleMiss_(false)
        {}


    // Member Functions

        //- Is there a hit
        __HOST____DEVICE__
        bool hit() const
        {
            return hit_;
        }

        //- Return hit point
        __HOST____DEVICE__
        const Point& hitPoint() const
        {/*
            if (!hit_)
            {
                FatalErrorIn("const Point& PointHit::hitPoint() const")
                    << "requested a hit point for a miss"
                    << abort(FatalError);
            }
*/
            return hitPoint_;
        }

        //- Return distance to hit
        __HOST____DEVICE__
        scalar distance() const
        {
            return distance_;
        }

        //- Return miss point
        __HOST____DEVICE__
        const Point& missPoint() const
        {/*
            if (hit_)
            {
                FatalErrorIn("const Point& PointHit::missPoint() const")
                    << "requested a miss point for a hit"
                    << abort(FatalError);
            }
*/
            return hitPoint_;
        }

        //- Return point with no checking
        __HOST____DEVICE__
        const Point& rawPoint() const
        {
            return hitPoint_;
        }

        //- Is this an eligible miss
        __HOST____DEVICE__
        bool eligibleMiss() const
        {
            return eligibleMiss_;
        }

        __HOST____DEVICE__
        void setHit()
        {
            hit_ = true;
            eligibleMiss_ = false;
        }

        __HOST____DEVICE__
        void setMiss(const bool eligible)
        {
            hit_ = false;
            eligibleMiss_ = eligible;
        }

        __HOST____DEVICE__
        void setPoint(const Point& p)
        {
            hitPoint_ = p;
        }

        __HOST____DEVICE__
        void setDistance(const scalar d)
        {
            distance_ = d;
        }


    // Ostream operator

        friend Ostream& operator<< <Point>
        (
            Ostream& os,
            const PointHit<Point>& b
        );
};


template<class Point>
inline Ostream& operator<<(Ostream& os, const PointHit<Point>& b)
{
    os  << b.hit() << token::SPACE
        << b.rawPoint() << token::SPACE
        << b.distance() << token::SPACE
        << b.eligibleMiss();

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
