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

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDistReflection.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);

    struct IDDESDeltaCalcDeltaFunctor
    {
        const label* cellFaces;
        const vector* faceCentres;
  
        IDDESDeltaCalcDeltaFunctor
        (
            const label* _cellFaces,
            const vector* _faceCentres
        ):
            cellFaces(_cellFaces),
            faceCentres(_faceCentres)
        {}

        __HOST____DEVICE__
        scalar operator()(const cellData& cell, const vector& nCell)
        {
            scalar deltaMaxTmp = 0.0;
            
            label start = cell.getStart();
            label end = start + cell.nFaces();

            for(label i = start; i < end; i++)
            {
                label faceI = cellFaces[i];
                const point& faceCentreI = faceCentres[faceI];

                for(label j = start; j < end; j++)
                {
                    label faceJ = cellFaces[j];
                    const point& faceCentreJ = faceCentres[faceJ];
                    scalar tmp = (faceCentreJ - faceCentreI) & nCell;
                    if (tmp > deltaMaxTmp)
                    {
                        deltaMaxTmp = tmp;
                    }
                }
            }

            return deltaMaxTmp;
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IDDESDelta::calcDelta()
{
    const volScalarField& hmax = hmax_();

    // initialise wallNorm
    wallDistReflection wallNorm(mesh());

    const volVectorField& n = wallNorm.n();

    tmp<volScalarField> tfaceToFacenMax
    (
        new volScalarField
        (
            IOobject
            (
                "faceToFaceMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    scalargpuField& faceToFacenMax = tfaceToFacenMax().internalField();

    const cellDatagpuList& cells = mesh().getCells();
    const labelgpuList& cellFaces = mesh().getCellFaces();

    const vectorgpuField& faceCentres = mesh().getFaceCentres();

    thrust::transform
    (
        cells.begin(),
        cells.end(),
        n.getField().begin(),
        faceToFacenMax.begin(),
        IDDESDeltaCalcDeltaFunctor
        (
            cellFaces.data(),
            faceCentres.data()
        )
    );

    label nD = mesh().nGeometricD();

    if (nD == 2)
    {
        WarningIn("IDDESDelta::calcDelta()")
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorIn("IDDESDelta::calcDelta()")
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    delta_.internalField() =
      ( 
        deltaCoeff_
       *min
        (
            max
            (
                max
                (
                    cw_*wallDist(mesh()).y(),
                    cw_*hmax
                ),
                tfaceToFacenMax
            ),
            hmax
        )
      )().getField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IDDESDelta::IDDESDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    hmax_(LESdelta::New("hmax", mesh, dd.parent())),
    deltaCoeff_(readScalar(dd.subDict(type()+"Coeffs").lookup("deltaCoeff"))),
    cw_(0.15)
{
    dd.subDict(type() + "Coeffs").readIfPresent("cw", cw_);
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IDDESDelta::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void Foam::IDDESDelta::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
        hmax_().correct();
    }
}


// ************************************************************************* //
