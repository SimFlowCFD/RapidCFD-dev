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

#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(maxDeltaxyz, 0);
addToRunTimeSelectionTable(LESdelta, maxDeltaxyz, dictionary);

struct maxDeltaxyzCalcDeltaFunctor
{
    const scalar deltaCoeff_;
    const label* cellFaces;
    const point* faceC;
    const vector* faceN;

    maxDeltaxyzCalcDeltaFunctor
    (
        const scalar deltaCoeff,
        const label* _cellFaces,
        const point* _faceC,
        const vector* _faceN
    ):
        deltaCoeff_(deltaCoeff),
        cellFaces(_cellFaces),
        faceC(_faceC),
        faceN(_faceN)
    {}

    __HOST____DEVICE__
    scalar operator()(const cellData& cell, const point& cc)
    {
        scalar deltaMaxTmp = 0.0;

        label start = cell.getStart();
        label end = start + cell.nFaces();

        for(label i = start; i < end; i++)
        {
            label faceI = cellFaces[i];

            const point& fc = faceC[faceI];
            const vector& n = faceN[faceI];

            scalar tmp = magSqr(n*(n & (fc - cc)));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }

        }

        return deltaCoeff_*sqrt(deltaMaxTmp);
    }
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void maxDeltaxyz::calcDelta()
{
    label nD = mesh().nGeometricD();

    tmp<volScalarField> hmax
    (
        new volScalarField
        (
            IOobject
            (
                "hmax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellDatagpuList& cells = mesh().getCells();
    const labelgpuList& cellFaces = mesh().getCellFaces();

    const vectorgpuField& cellC = mesh().getCellCentres();
    const vectorgpuField& faceC = mesh().getFaceCentres();
    const vectorgpuField faceN(mesh().getFaceAreas()/mag(mesh().getFaceAreas()));

    thrust::transform
    (
        cells.begin(),
        cells.end(),
        cellC.begin(),
        hmax().getField().begin(),
        maxDeltaxyzCalcDeltaFunctor
        (
            deltaCoeff_,
            cellFaces.data(),
            faceC.data(),
            faceN.data()
        )
    );

    if (nD == 3)
    {
        delta_.internalField() = hmax();
    }
    else if (nD == 2)
    {
        WarningIn("maxDeltaxyz::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        delta_.internalField() = hmax();
    }
    else
    {
        FatalErrorIn("maxDeltaxyz::calcDelta()")
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxDeltaxyz::maxDeltaxyz
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    deltaCoeff_(readScalar(dd.subDict(type() + "Coeffs").lookup("deltaCoeff")))
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void maxDeltaxyz::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void maxDeltaxyz::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
