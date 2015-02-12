/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "hExponentialThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hExponentialThermo<EquationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0.0)
    {
        FatalErrorIn
        (
            "hExponentialThermo<EquationOfState>::checkT(const scalar T) const"
        )   << "attempt to use hExponentialThermo<EquationOfState>"
               " out of temperature range "
            << T
            << abort(FatalError);
    }
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::integrateCp
(
    const scalar T
) const
{
    return
    (
        c0_*pow(T, n0_ + 1.0)/(pow(Tref_, n0_)*(n0_ + 1.0))
    );
}


template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState>::hExponentialThermo
(
    const word& name,
    const hExponentialThermo& jt
)
:
    EquationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    Hf_(jt.Hf_)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState>::hExponentialThermo
(
    const EquationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar Hf
)
:
    EquationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    Hf_(Hf)
{}

template<class EquationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<EquationOfState> >
Foam::hExponentialThermo<EquationOfState>::clone() const
{
    return autoPtr<hExponentialThermo<EquationOfState> >
    (
        new hExponentialThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<EquationOfState> >
Foam::hExponentialThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<hExponentialThermo<EquationOfState> >
    (
        new hExponentialThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hExponentialThermo<EquationOfState> >
Foam::hExponentialThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hExponentialThermo<EquationOfState> >
    (
        new hExponentialThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::cp
(
    const scalar p, const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_);
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    scalar hOffset = integrateCp(298.15);

    return
    (
        (integrateCp(T) + Hf_ - hOffset)
    );
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    scalar hOffset = integrateCp(298.15);
    return (integrateCp(T) - hOffset);
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::hc() const
{
    return Hf_;
}


template<class EquationOfState>
__HOST____DEVICE__
inline Foam::scalar Foam::hExponentialThermo<EquationOfState>::s
(
    const scalar p, const scalar T
) const
{
    #ifndef __CUDACC__
    notImplemented
    (
        "scalar hExponentialThermo<EquationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    #endif
    return T;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hExponentialThermo<EquationOfState>::operator+=
(
    const hExponentialThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(ct);
    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
    c0_ = molr1*c0_ + molr2*ct.c0_;
    n0_ = molr1*n0_ + molr2*ct.n0_;
    Tref_ = molr1*Tref_ + molr2*ct.Tref_;
}


template<class EquationOfState>
inline void Foam::hExponentialThermo<EquationOfState>::operator-=
(
    const hExponentialThermo<EquationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
    c0_ = (molr1*c0_ - molr2*ct.c0_);
    n0_ = (molr1*n0_ - molr2*ct.n0_);
    Tref_ = (molr1*Tref_ - molr2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState> Foam::operator+
(
    const hExponentialThermo<EquationOfState>& ct1,
    const hExponentialThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return hExponentialThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.c0_
      + ct2.nMoles()/eofs.nMoles()*ct2.c0_,
        ct1.nMoles()/eofs.nMoles()*ct1.n0_
      + ct2.nMoles()/eofs.nMoles()*ct2.n0_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tref_
      + ct2.nMoles()/eofs.nMoles()*ct2.Tref_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState> Foam::operator-
(
    const hExponentialThermo<EquationOfState>& ct1,
    const hExponentialThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    return hExponentialThermo<EquationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.c0_
      - ct2.nMoles()/eofs.nMoles()*ct2.c0_,
        ct1.nMoles()/eofs.nMoles()*ct1.n0_
      - ct2.nMoles()/eofs.nMoles()*ct2.n0_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tref_
      - ct2.nMoles()/eofs.nMoles()*ct2.Tref_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hExponentialThermo<EquationOfState>& ct
)
{
    return hExponentialThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.Hf_
    );
}


template<class EquationOfState>
inline Foam::hExponentialThermo<EquationOfState> Foam::operator==
(
    const hExponentialThermo<EquationOfState>& ct1,
    const hExponentialThermo<EquationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
