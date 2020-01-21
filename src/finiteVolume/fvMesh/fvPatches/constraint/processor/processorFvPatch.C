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

#include "processorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, processorFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

namespace Foam
{
struct processorFvPatchMakeWeights
{
    template<class Tuple>
    __host__ __device__
    scalar operator()(const Tuple& t)
    {
        const vector& nfa = thrust::get<0>(t);
        const vector& nfc = thrust::get<1>(t);
        const vector& nfcc = thrust::get<2>(t);
        const vector& nf = thrust::get<3>(t);
        const vector& delta = thrust::get<4>(t);

        // The face normals point in the opposite direction on the other side
        vector nfaHat = nfa/(mag(nfa)+VSMALL);
        scalar nfccn = (nfaHat  & (nfc - nfcc));

        return nfccn/((nf&delta)+nfccn);
    }
};
}

void Foam::processorFvPatch::makeWeights(scalargpuField& w) const
{
    if (Pstream::parRun())
    {
        tmp<vectorgpuField> tdelt = coupledFvPatch::delta();
        const vectorgpuField& delt = tdelt();

        tmp<vectorgpuField> tnf = nf();
        const vectorgpuField& nf = tnf();

        auto iter = thrust::make_zip_iterator(thrust::make_tuple(
            procPolyPatch_.getNeighbFaceAreas().begin(),
            procPolyPatch_.getNeighbFaceCentres().begin(),
            procPolyPatch_.getNeighbFaceCellCentres().begin(),
            nf.begin(),
            delt.begin()
        ));
        thrust::transform(iter, iter+size(), w.begin(),
             processorFvPatchMakeWeights()
        );
       /* scalargpuField neighbFaceCentresCn
        (
            (
                procPolyPatch_.getNeighbFaceAreas()
               /(mag(procPolyPatch_.getNeighbFaceAreas()) + VSMALL)
            )
          & (
              procPolyPatch_.getNeighbFaceCentres()
            - procPolyPatch_.getNeighbFaceCellCentres())
        );

        w = neighbFaceCentresCn
           /((nf()&coupledFvPatch::delta()) + neighbFaceCentresCn);*/
    }
    else
    {
        w = 1.0;
    }
}


Foam::tmp<Foam::vectorgpuField> Foam::processorFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                coupledFvPatch::delta()
              - (
                    procPolyPatch_.getNeighbFaceCentres()
                  - procPolyPatch_.getNeighbFaceCellCentres()
                );
        }
        else
        {
            return
                coupledFvPatch::delta()
              - transform
                (
                    getForwardT(),
                    (
                        procPolyPatch_.getNeighbFaceCentres()
                      - procPolyPatch_.getNeighbFaceCellCentres()
                    )
                );
        }
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::processorFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::processorFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    tmp<labelField> tmpf(new labelField(this->size()));
    receive(commsType, tmpf());
    return tmpf;
}


// ************************************************************************* //
