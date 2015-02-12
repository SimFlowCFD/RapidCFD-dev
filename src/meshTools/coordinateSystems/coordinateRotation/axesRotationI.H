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

\*---------------------------------------------------------------------------*/

#include "axesRotation.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

__HOST____DEVICE__
inline Foam::vector Foam::axesRotation::transform(const vector& st) const
{
    return (R_ & st);
}

__HOST____DEVICE__
inline Foam::vector Foam::axesRotation::invTransform(const vector& st) const
{
    return (Rtr_ & st);
}

__HOST____DEVICE__
inline Foam::tensor Foam::axesRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}

__HOST____DEVICE__
inline Foam::symmTensor Foam::axesRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}

// ************************************************************************* //
