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

#include "tensor.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* const tensor::typeName = "tensor";

    template<>
    const char* tensor::componentNames[] =
    {
        "xx", "xy", "xz",
        "yx", "yy", "yz",
        "zx", "zy", "zz"
    };

    template<>
    const tensor tensor::zero
    (
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    );

    template<>
    const tensor tensor::one
    (
        1, 1, 1,
        1, 1, 1,
        1, 1, 1
    );

    template<>
    const tensor tensor::max
    (
        VGREAT, VGREAT, VGREAT,
        VGREAT, VGREAT, VGREAT,
        VGREAT, VGREAT, VGREAT
    );

    template<>
    const tensor tensor::min
    (
        -VGREAT, -VGREAT, -VGREAT,
        -VGREAT, -VGREAT, -VGREAT,
        -VGREAT, -VGREAT, -VGREAT
    );

    template<>
    const tensor tensor::I
    (
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    );
}

// ************************************************************************* //
