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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

__HOST____DEVICE__
inline Foam::quaternion::quaternion()
{}

__HOST____DEVICE__
inline Foam::quaternion::quaternion(const scalar w, const vector& v)
:
    w_(w),
    v_(v)
{}

__HOST____DEVICE__
inline Foam::quaternion::quaternion(const vector& d, const scalar theta)
:
    w_(cos(0.5*theta)),
    v_((sin(0.5*theta)/mag(d))*d)
{}

__HOST____DEVICE__
inline Foam::quaternion::quaternion
(
    const vector& d,
    const scalar cosTheta,
    const bool normalized
)
{
    scalar cosHalfTheta2 = 0.5*(cosTheta + 1);
    w_ = sqrt(cosHalfTheta2);

    if (normalized)
    {
        v_ = sqrt(1 - cosHalfTheta2)*d;
    }
    else
    {
        v_ = (sqrt(1 - cosHalfTheta2)/mag(d))*d;
    }
}

__HOST____DEVICE__
inline Foam::quaternion::quaternion(const scalar w)
:
    w_(w),
    v_(0,0,0)
{}

__HOST____DEVICE__
inline Foam::quaternion::quaternion(const vector& v)
:
    w_(0),
    v_(v)
{}

__HOST____DEVICE__
inline Foam::quaternion::quaternion
(
    const scalar angleX,
    const scalar angleY,
    const scalar angleZ
)
{
    operator=(quaternion(vector(1, 0, 0), angleX));
    operator*=(quaternion(vector(0, 1, 0), angleY));
    operator*=(quaternion(vector(0, 0, 1), angleZ));
}

__HOST____DEVICE__
inline Foam::quaternion::quaternion
(
    const tensor& rotationTensor
)
{
    scalar trace =
        rotationTensor.xx()
      + rotationTensor.yy()
      + rotationTensor.zz();

    if (trace > 0)
    {
        scalar s = 0.5/Foam::sqrt(trace + 1.0);

        w_ = 0.25/s;
        v_[0] = (rotationTensor.zy() - rotationTensor.yz())*s;
        v_[1] = (rotationTensor.xz() - rotationTensor.zx())*s;
        v_[2] = (rotationTensor.yx() - rotationTensor.xy())*s;
    }
    else
    {
        if
        (
            rotationTensor.xx() > rotationTensor.yy()
         && rotationTensor.xx() > rotationTensor.zz()
        )
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.xx()
              - rotationTensor.yy()
              - rotationTensor.zz()
            );

            w_ = (rotationTensor.zy() - rotationTensor.yz())/s;
            v_[0] = 0.25*s;
            v_[1] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v_[2] = (rotationTensor.xz() + rotationTensor.zx())/s;
        }
        else if
        (
            rotationTensor.yy() > rotationTensor.zz()
        )
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.yy()
              - rotationTensor.xx()
              - rotationTensor.zz()
            );

            w_ = (rotationTensor.xz() - rotationTensor.xz())/s;
            v_[0] = (rotationTensor.xy() + rotationTensor.yx())/s;
            v_[1] = 0.25*s;
            v_[2] = (rotationTensor.yz() + rotationTensor.zy())/s;
        }
        else
        {
            scalar s = 2.0*Foam::sqrt
            (
                1.0
              + rotationTensor.zz()
              - rotationTensor.xx()
              - rotationTensor.yy()
            );

            w_ = (rotationTensor.yx() - rotationTensor.xy())/s;
            v_[0] = (rotationTensor.xz() + rotationTensor.zx())/s;
            v_[1] = (rotationTensor.yz() + rotationTensor.zy())/s;
            v_[2] = 0.25*s;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline Foam::scalar Foam::quaternion::w() const
{
    return w_;
}


__HOST____DEVICE__
inline const Foam::vector& Foam::quaternion::v() const
{
    return v_;
}


__HOST____DEVICE__
inline Foam::scalar& Foam::quaternion::w()
{
    return w_;
}


__HOST____DEVICE__
inline Foam::vector& Foam::quaternion::v()
{
    return v_;
}


__HOST____DEVICE__
inline Foam::quaternion Foam::quaternion::normalized() const
{
    return *this/mag(*this);
}


__HOST____DEVICE__
inline void Foam::quaternion::normalize()
{
    operator/=(mag(*this));
}


__HOST____DEVICE__
inline Foam::quaternion Foam::quaternion::mulq0v(const vector& u) const
{
    return quaternion(-(v() & u), w()*u + (v() ^ u));
}


__HOST____DEVICE__
inline Foam::vector Foam::quaternion::transform(const vector& u) const
{
    return (mulq0v(u)*conjugate(*this)).v();
}


__HOST____DEVICE__
inline Foam::vector Foam::quaternion::invTransform(const vector& u) const
{
    return (conjugate(*this).mulq0v(u)*(*this)).v();
}


__HOST____DEVICE__
inline Foam::quaternion Foam::quaternion::transform(const quaternion& q) const
{
    return Foam::normalize((*this)*q);
}


__HOST____DEVICE__
inline Foam::quaternion Foam::quaternion::invTransform
(
    const quaternion& q
) const
{
    return Foam::normalize(conjugate(*this)*q);
}


__HOST____DEVICE__
inline Foam::tensor Foam::quaternion::R() const
{
    scalar w2 = sqr(w());
    scalar x2 = sqr(v().x());
    scalar y2 = sqr(v().y());
    scalar z2 = sqr(v().z());

    scalar txy = 2*v().x()*v().y();
    scalar twz = 2*w()*v().z();
    scalar txz = 2*v().x()*v().z();
    scalar twy = 2*w()*v().y();
    scalar tyz = 2*v().y()*v().z();
    scalar twx = 2*w()*v().x();

    return tensor
    (
        w2 + x2 - y2 - z2,  txy - twz,          txz + twy,
        txy + twz,          w2 - x2 + y2 - z2,  tyz - twx,
        txz - twy,          tyz + twx,          w2 - x2 - y2 + z2
    );
}


__HOST____DEVICE__
inline Foam::vector Foam::quaternion::eulerAngles(const quaternion& q) const
{
    vector angles(vector::zero);

    const scalar& w = q.w();
    const vector& v = q.v();

    angles[0] = Foam::atan2
    (
        2*(w*v.x() + v.y()*v.z()),
        1 - 2*(sqr(v.x()) + sqr(v.y()))
    );
    angles[1] = Foam::asin(2*(w*v.y() - v.z()*v.x()));
    angles[2] = Foam::atan2
    (
        2*(w*v.z() + v.x()*v.y()),
        1 - 2*(sqr(v.y()) + sqr(v.z()))
    );

    // Convert to degrees
//    forAll(angles, aI)
//    {
//        angles[aI] = Foam::radToDeg(angles[aI]);
//    }

    return angles;
}


__HOST____DEVICE__
inline Foam::quaternion Foam::slerp
(
    const quaternion& qa,
    const quaternion& qb,
    const scalar t
)
{
    label sign = 1;

    if ((qa & qb) < 0)
    {
        sign = -1;
    }

    return qa*pow((inv(qa)*sign*qb), t);
}

__HOST____DEVICE__
inline Foam::quaternion Foam::exp(const quaternion& q)
{
    const scalar magV = mag(q.v());

    if (magV == 0)
    {
        return quaternion(1, vector(0,0,0));
    }

    const scalar expW = exp(q.w());

    return quaternion
    (
        expW*cos(magV),
        expW*sin(magV)*q.v()/magV
    );
}

__HOST____DEVICE__
inline Foam::quaternion Foam::pow(const quaternion& q, const label power)
{
    const scalar magQ = mag(q);
    const scalar magV = mag(q.v());

    quaternion powq(q.v());

    if (magV != 0 && magQ != 0)
    {
        powq /= magV;
        powq *= power*acos(q.w()/magQ);
    }

    return pow(magQ, power)*exp(powq);
}

__HOST____DEVICE__
inline Foam::quaternion Foam::pow(const quaternion& q, const scalar power)
{
    const scalar magQ = mag(q);
    const scalar magV = mag(q.v());

    quaternion powq(q.v());

    if (magV != 0 && magQ != 0)
    {
        powq /= magV;
        powq *= power*acos(q.w()/magQ);
    }

    return pow(magQ, power)*exp(powq);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline void Foam::quaternion::operator=(const quaternion& q)
{
    w_ = q.w_;
    v_ = q.v_;
}

__HOST____DEVICE__
inline void Foam::quaternion::operator+=(const quaternion& q)
{
    w_ += q.w_;
    v_ += q.v_;
}

__HOST____DEVICE__
inline void Foam::quaternion::operator-=(const quaternion& q)
{
    w_ -= q.w_;
    v_ -= q.v_;
}

__HOST____DEVICE__
inline void Foam::quaternion::operator*=(const quaternion& q)
{
    scalar w0 = w();
    w() = w()*q.w() - (v() & q.v());
    v() = w0*q.v() + q.w()*v() + (v() ^ q.v());
}

inline void Foam::quaternion::operator/=(const quaternion& q)
{
    return operator*=(inv(q));
}


__HOST____DEVICE__
inline void Foam::quaternion::operator=(const scalar s)
{
    w_ = s;
}


__HOST____DEVICE__
inline void Foam::quaternion::operator=(const vector& v)
{
    v_ = v;
}


__HOST____DEVICE__
inline void Foam::quaternion::operator*=(const scalar s)
{
    w_ *= s;
    v_ *= s;
}

__HOST____DEVICE__
inline void Foam::quaternion::operator/=(const scalar s)
{
    w_ /= s;
    v_ /= s;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline Foam::scalar Foam::magSqr(const quaternion& q)
{
    return magSqr(q.w()) + magSqr(q.v());
}


__HOST____DEVICE__
inline Foam::scalar Foam::mag(const quaternion& q)
{
    return sqrt(magSqr(q));
}


__HOST____DEVICE__
inline Foam::quaternion Foam::conjugate(const quaternion& q)
{
    return quaternion(q.w(), -q.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::inv(const quaternion& q)
{
    scalar magSqrq = magSqr(q);
    return quaternion(q.w()/magSqrq, -q.v()/magSqrq);
}


__HOST____DEVICE__
inline Foam::quaternion Foam::normalize(const quaternion& q)
{
    return q/mag(q);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

__HOST____DEVICE__
inline bool Foam::operator==(const quaternion& q1, const quaternion& q2)
{
    return (equal(q1.w(), q2.w()) && equal(q1.v(), q2.v()));
}


__HOST____DEVICE__
inline bool Foam::operator!=(const quaternion& q1, const quaternion& q2)
{
    return !operator==(q1, q2);
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator+
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion(q1.w() + q2.w(), q1.v() + q2.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator-(const quaternion& q)
{
    return quaternion(-q.w(), -q.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator-
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion(q1.w() - q2.w(), q1.v() - q2.v());
}


__HOST____DEVICE__
inline Foam::scalar Foam::operator&(const quaternion& q1, const quaternion& q2)
{
    return q1.w()*q2.w() + (q1.v() & q2.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator*
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return quaternion
    (
        q1.w()*q2.w() - (q1.v() & q2.v()),
        q1.w()*q2.v() + q2.w()*q1.v() + (q1.v() ^ q2.v())
    );
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator/
(
    const quaternion& q1,
    const quaternion& q2
)
{
    return q1*inv(q2);
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator*(const scalar s, const quaternion& q)
{
    return quaternion(s*q.w(), s*q.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator*(const quaternion& q, const scalar s)
{
    return quaternion(s*q.w(), s*q.v());
}


__HOST____DEVICE__
inline Foam::quaternion Foam::operator/(const quaternion& q, const scalar s)
{
    return quaternion(q.w()/s, q.v()/s);
}


// ************************************************************************* //
