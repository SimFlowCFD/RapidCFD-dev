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

#include "movingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "faceFunctors.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallVelocityFvPatchVectorField::
movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_("U")
{}


Foam::movingWallVelocityFvPatchVectorField::
movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


Foam::movingWallVelocityFvPatchVectorField::
movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    fvPatchVectorField::operator=(vectorgpuField("value", dict, p.size()));
}


Foam::movingWallVelocityFvPatchVectorField::
movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf),
    UName_(mwvpvf.UName_)
{}


Foam::movingWallVelocityFvPatchVectorField::
movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF),
    UName_(mwvpvf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = dimensionedInternalField().mesh();

    if (mesh.moving())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointgpuField& oldPoints = mesh.getOldPoints();

        const faceDatagpuList& faces = pp.getFaces();
        const labelgpuList& faceNodes = pp.getFaceNodes();

        vectorgpuField oldFc(pp.size());

        thrust::transform
        (
            faces.begin(),
            faces.end(),
            oldFc.begin(),
            faceCentreFunctor
            (
                faceNodes.data(),
                oldPoints.data()
            )
        );

        const scalar deltaT = mesh.time().deltaTValue();

        const vectorgpuField Up((pp.getFaceCentres() - oldFc)/deltaT);

        const volVectorField& U = db().lookupObject<volVectorField>(UName_);
        scalargpuField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );

        const vectorgpuField n(p.nf());
        const scalargpuField& magSf = p.magSf();
        tmp<scalargpuField> Un = phip/(magSf + VSMALL);


        vectorgpuField::operator=(Up + n*(Un - (n & Up)));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
