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

Typedef
    Foam::scalar

Description
    Single floating point number identical to float or double depending on
    whether WM_SP or WM_DP is defined.

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "floatScalar.H"
#include "doubleScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;

    const scalar GREAT = floatScalarGREAT;
    const scalar VGREAT = floatScalarVGREAT;
    const scalar ROOTVGREAT = floatScalarROOTVGREAT;
    const scalar SMALL = floatScalarSMALL;
    const scalar VSMALL = floatScalarVSMALL;
    const scalar ROOTVSMALL = floatScalarROOTVSMALL;

    __constant__ const scalar GPU_GREAT = floatScalarGREAT;
    __constant__ const scalar GPU_VGREAT = floatScalarVGREAT;
    __constant__ const scalar GPU_ROOTVGREAT = floatScalarROOTVGREAT;
    __constant__ const scalar GPU_SMALL = floatScalarSMALL;
    __constant__ const scalar GPU_VSMALL = floatScalarVSMALL;
    __constant__ const scalar GPU_ROOTVSMALL = floatScalarROOTVSMALL;

    scalar readScalar(Istream& is);
}

#elif defined(WM_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;

    const scalar GREAT = doubleScalarGREAT;
    const scalar VGREAT = doubleScalarVGREAT;
    const scalar ROOTVGREAT = doubleScalarROOTVGREAT;
    const scalar SMALL = doubleScalarSMALL;
    const scalar VSMALL = doubleScalarVSMALL;
    const scalar ROOTVSMALL = doubleScalarROOTVSMALL;

    __constant__ const scalar GPU_GREAT = doubleScalarGREAT;
    __constant__ const scalar GPU_VGREAT = doubleScalarVGREAT;
    __constant__ const scalar GPU_ROOTVGREAT = doubleScalarROOTVGREAT;
    __constant__ const scalar GPU_SMALL = doubleScalarSMALL;
    __constant__ const scalar GPU_VSMALL = doubleScalarVSMALL;
    __constant__ const scalar GPU_ROOTVSMALL = doubleScalarROOTVSMALL;

    scalar readScalar(Istream& is);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
