/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "cyclicACMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
    
    struct cyclicACMIFvPatchMakeWeightsFunctor
    {
        __HOST____DEVICE__
        scalar operator () (const scalar& deltas, const scalar& nbrDeltas)
        {
            return nbrDeltas/(deltas + nbrDeltas);
        }
    };
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicACMIFvPatch::updateAreas() const
{
    if (cyclicACMIPolyPatch_.updated())
    {
        // Set Sf and magSf for both sides' coupled and non-overlapping patches

        // owner couple
        const_cast<vectorgpuField&>(Sf()) = patch().faceAreas();
        const_cast<scalargpuField&>(magSf()) = mag(patch().faceAreas());

        // owner non-overlapping
        const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
        const_cast<vectorgpuField&>(nonOverlapPatch.Sf()) =
            nonOverlapPatch.patch().faceAreas();
        const_cast<scalargpuField&>(nonOverlapPatch.magSf()) =
            mag(nonOverlapPatch.patch().faceAreas());

        // neighbour couple
        const cyclicACMIFvPatch& nbrACMI = neighbPatch();
        const_cast<vectorgpuField&>(nbrACMI.Sf()) =
            nbrACMI.patch().faceAreas();
        const_cast<scalargpuField&>(nbrACMI.magSf()) =
            mag(nbrACMI.patch().faceAreas());

        // neighbour non-overlapping
        const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();
        const_cast<vectorgpuField&>(nbrNonOverlapPatch.Sf()) =
            nbrNonOverlapPatch.patch().faceAreas();
        const_cast<scalargpuField&>(nbrNonOverlapPatch.magSf()) =
            mag(nbrNonOverlapPatch.patch().faceAreas());

        // set the updated flag
        cyclicACMIPolyPatch_.setUpdated(false);
    }
}


void Foam::cyclicACMIFvPatch::makeWeights(scalargpuField& w) const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        const fvPatch& nbrPatchNonOverlap = nonOverlapPatch();

        const scalargpuField deltas(nf() & coupledFvPatch::delta());

        const scalargpuField nbrDeltas
        (
            interpolate
            (
                nbrPatch.nf() & nbrPatch.coupledFvPatch::delta(),
                nbrPatchNonOverlap.nf() & nbrPatchNonOverlap.delta()
            )
        );

        thrust::transform
        (
            deltas.begin(),
            deltas.end(),
            nbrDeltas.begin(),
            w.begin(),
            cyclicACMIFvPatchMakeWeightsFunctor()
        );
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicACMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


Foam::tmp<Foam::vectorgpuField> Foam::cyclicACMIFvPatch::delta() const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatchCoupled = neighbFvPatch();
        const fvPatch& nbrPatchNonOverlap = nonOverlapPatch();

        const vectorgpuField patchD(coupledFvPatch::delta());

        vectorgpuField nbrPatchD
        (
            interpolate
            (
                nbrPatchCoupled.coupledFvPatch::delta(),
                nbrPatchNonOverlap.delta()
            )
        );

        const vectorgpuField nbrPatchD0
        (
            interpolate
            (
                vectorgpuField(nbrPatchCoupled.size(), vector::zero),
                nbrPatchNonOverlap.delta()()
            )
        );

        nbrPatchD -= nbrPatchD0;

        tmp<vectorgpuField> tpdv(new vectorgpuField(patchD.size()));
        vectorgpuField& pdv = tpdv();

        // do the transformation if necessary
        if (parallel())
        {
            thrust::transform
            (
                patchD.begin(),
                patchD.end(),
                nbrPatchD.begin(),
                pdv.begin(),
                subtractOperatorFunctor<vector,vector,vector>()
            );
        }
        else
        {
            tensor t = forwardT()[0];
			
            thrust::transform
            (
                patchD.begin(),
                patchD.end(),
                thrust::make_transform_iterator
                (
                    nbrPatchD.begin(),
                    transformBinaryFunctionSFFunctor<tensor,vector,vector>(t)
                ),
                pdv.begin(),
                subtractOperatorFunctor<vector,vector,vector>()
            );
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
