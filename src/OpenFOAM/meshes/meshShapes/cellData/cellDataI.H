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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
__HOST____DEVICE__
inline Foam::cellData::cellData()
{}


// Construct given size
__HOST____DEVICE__
inline Foam::cellData::cellData(label st, label si)
:
    start(st),
    size(si)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
__HOST____DEVICE__
inline bool Foam::operator==(const cellData& a, const cellData& b)
{
    return (a.start == b.start);
}

__HOST____DEVICE__
inline bool Foam::operator!=(const cellData& a, const cellData& b)
{
    return (!(a == b));
}

// ************************************************************************* //
