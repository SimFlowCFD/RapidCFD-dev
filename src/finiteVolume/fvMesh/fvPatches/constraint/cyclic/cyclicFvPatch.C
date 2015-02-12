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

#include "cyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicFvPatch, polyPatch);

    struct cyclicMakeWeightsFunctor : public std::binary_function<scalar,scalar,scalar>{
         __HOST____DEVICE__
         scalar operator()(const scalar& di,const scalar& dni){
             return dni/(di + dni);
         }
    };
    
    struct cyclicDeltaFunctor : public std::binary_function<vector,vector,vector>{
		const tensor t;
		cyclicDeltaFunctor(const tensor _t):t(_t){}
		__HOST____DEVICE__
		vector operator()(const vector& ddi, const vector& dni){
			return ddi - transform(t, dni);
		}
	};
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicFvPatch::makeWeights(scalargpuField& w) const
{
    const cyclicFvPatch& nbrPatch = neighbFvPatch();

    const scalargpuField deltas(nf()&coupledFvPatch::delta());
    const scalargpuField nbrDeltas(nbrPatch.nf()&nbrPatch.coupledFvPatch::delta());
/*
    forAll(deltas, facei)
    {
        scalar di = deltas[facei];
        scalar dni = nbrDeltas[facei];

        w[facei] = dni/(di + dni);
    }
*/
    thrust::transform(deltas.begin(),deltas.end(),nbrDeltas.begin(),w.begin(),
                      cyclicMakeWeightsFunctor());
}


Foam::tmp<Foam::vectorgpuField> Foam::cyclicFvPatch::delta() const
{
    const vectorgpuField patchD(coupledFvPatch::delta());
    const vectorgpuField nbrPatchD(neighbFvPatch().coupledFvPatch::delta());

    tmp<vectorgpuField> tpdv(new vectorgpuField(patchD.size()));
    vectorgpuField& pdv = tpdv();

    // To the transformation if necessary
    if (parallel())
    {
/*
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - dni;
        }
*/
        thrust::transform(patchD.begin(),patchD.end(),nbrPatchD.begin(),pdv.begin(),
                          subtractOperatorFunctor<vector,vector,vector>());
    }
    else
    {
/*
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - transform(forwardT()[0], dni);
        }
*/
        thrust::transform(patchD.begin(),patchD.end(),nbrPatchD.begin(),pdv.begin(),
                          cyclicDeltaFunctor(forwardT()[0]));
    }

    return tpdv;
}


Foam::tmp<Foam::labelField> Foam::cyclicFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
