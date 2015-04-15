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

#include "cyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicAMIFvPatch, polyPatch);
    
    struct cyclicAMIFvPatchMakeWeightsFunctor
    {
        __HOST____DEVICE__
        scalar operator () (const scalar& deltas, const scalar& nbrDeltas)
        {
            return nbrDeltas/(deltas + nbrDeltas);
        }
    };
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicAMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


void Foam::cyclicAMIFvPatch::makeWeights(scalargpuField& w) const
{
    if (coupled())
    {
        const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

        const scalargpuField deltas(nf() & coupledFvPatch::delta());

        tmp<scalargpuField> tnbrDeltas;
        if (applyLowWeightCorrection())
        {
            tnbrDeltas =
                interpolate
                (
                    nbrPatch.nf() & nbrPatch.coupledFvPatch::delta(),
                    scalargpuField(this->size(), 1.0)
                );
        }
        else
        {
            tnbrDeltas =
                interpolate(nbrPatch.nf() & nbrPatch.coupledFvPatch::delta());
        }

        const scalargpuField& nbrDeltas = tnbrDeltas();

        thrust::transform
        (
            deltas.begin(),
            deltas.end(),
            nbrDeltas.begin(),
            w.begin(),
            cyclicAMIFvPatchMakeWeightsFunctor()
        );
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


Foam::tmp<Foam::vectorgpuField> Foam::cyclicAMIFvPatch::delta() const
{
    const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

    if (coupled())
    {
        const vectorgpuField patchD(coupledFvPatch::delta());

        tmp<vectorgpuField> tnbrPatchD;
        if (applyLowWeightCorrection())
        {
            tnbrPatchD =
                interpolate
                (
                    nbrPatch.coupledFvPatch::delta(),
                    vectorgpuField(this->size(), vector::zero)
                );
        }
        else
        {
            tnbrPatchD = interpolate(nbrPatch.coupledFvPatch::delta());
        }

        const vectorgpuField& nbrPatchD = tnbrPatchD();

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


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
