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

#include "localAxesRotation.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "tensorIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localAxesRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        localAxesRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        localAxesRotation,
        objectRegistry
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::localAxesRotation::init
(
    const objectRegistry& obr,
    const List<label>& cells
)
{
    const polyMesh& mesh = refCast<const polyMesh>(obr);
    const vectorField& cc = mesh.cellCentres();

    if (cells.size())
    {
        Rptr_.reset(new tensorField(cells.size()));

        tensorField& R = Rptr_();
        forAll(cells, i)
        {
            label cellI = cells[i];
            vector dir = cc[cellI] - origin_;
            dir /= mag(dir) + VSMALL;

            R[i] = axesRotation(e3_, dir).R();
        }
    }
    else
    {
        Rptr_.reset(new tensorField(mesh.nCells()));

        tensorField& R = Rptr_();
        forAll(cc, cellI)
        {
            vector dir = cc[cellI] - origin_;
            dir /= mag(dir) + VSMALL;

            R[cellI] = axesRotation(e3_, dir).R();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localAxesRotation::localAxesRotation
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    Rptr_(),
    Rgpuptr_(),
    origin_(point::zero),
    e3_(vector::zero)
{
    // If origin is specified in the coordinateSystem
    if (dict.parent().found("origin"))
    {
        dict.parent().lookup("origin") >> origin_;
    }

    // rotation axis
    dict.lookup("e3") >> e3_;

    init(obr);
}


Foam::localAxesRotation::localAxesRotation
(
    const objectRegistry& obr,
    const vector& axis,
    const point& origin
)
:
    Rptr_(),
    Rgpuptr_(),
    origin_(origin),
    e3_(axis)
{
    init(obr);
}


Foam::localAxesRotation::localAxesRotation
(
    const objectRegistry& obr,
    const vector& axis,
    const point& origin,
    const List<label>& cells
)
:
    Rptr_(),
    Rgpuptr_(),
    origin_(origin),
    e3_(axis)
{
    init(obr, cells);
}


Foam::localAxesRotation::localAxesRotation(const dictionary& dict)
:
    Rptr_(),
    Rgpuptr_(),
    origin_(),
    e3_()
{
    FatalErrorIn("localAxesRotation(const dictionary&)")
        << " localAxesRotation can not be constructed from dictionary "
        << " use the construtctor : "
           "("
           "    const dictionary&, const objectRegistry&"
           ")"
        << exit(FatalIOError);
}


Foam::localAxesRotation::localAxesRotation(const tensorField& R)
:
    Rptr_(),
    Rgpuptr_(),
    origin_(vector::zero),
    e3_(vector::zero)
{
    Rptr_() = R;
    Rgpuptr_() = R;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::localAxesRotation::clear()
{
    if (!Rptr_.empty())
    {
        Rptr_.clear();
    }
    if (!Rgpuptr_.empty())
    {
        Rgpuptr_.clear();
    }
}


void Foam::localAxesRotation::updateCells
(
    const polyMesh& mesh,
    const labelgpuList& cells
)
{
    const vectorgpuField& cc = mesh.getCellCentres();
    tensorgpuField& R = Rgpuptr_();

    forAll(cells, i)
    {
        label cellI = cells[i];
        vector dir = cc[cellI] - origin_;
        dir /= mag(dir) + VSMALL;

        R[cellI] = axesRotation(e3_, dir).R();
    }
}


Foam::tmp<Foam::vectorField> Foam::localAxesRotation::transform
(
    const vectorField& vf
) const
{
    if (Rptr_->size() != vf.size())
    {
        FatalErrorIn
        (
            "tmp<vectorField> localAxesRotation::transform(const vectorField&)"
        )
            << "vectorField st has different size to tensorField "
            << abort(FatalError);
    }

    return (Rptr_() & vf);
}

__HOST____DEVICE__
Foam::vector Foam::localAxesRotation::transform(const vector& v) const
{
    #ifndef __CUDACC__
    notImplemented
    (
        "vector localAxesRotation::transform(const vector&) const"
    );
    return vector::zero;
    #else
    return vector(0,0,0);
    #endif
}


__HOST____DEVICE__
Foam::vector Foam::localAxesRotation::transform
(
    const vector& v,
    const label cmptI
) const
{
    return (Rptr_()[cmptI] & v);
}


Foam::tmp<Foam::vectorField> Foam::localAxesRotation::invTransform
(
    const vectorField& vf
) const
{
    return (Rptr_().T() & vf);
}


__HOST____DEVICE__
Foam::vector Foam::localAxesRotation::invTransform(const vector& v) const
{
    #ifndef __CUDACC__
    notImplemented
    (
        "vector localAxesRotation::invTransform(const vector&) const"
    );
    return vector::zero;
    #else
    return vector(0,0,0);
    #endif
}


__HOST____DEVICE__
Foam::vector Foam::localAxesRotation::invTransform
(
    const vector& v,
    const label cmptI
) const
{
    return (Rptr_()[cmptI].T() & v);
}


Foam::tmp<Foam::tensorField> Foam::localAxesRotation::transformTensor
(
    const tensorField& tf
) const
{
    if (Rptr_->size() != tf.size())
    {
        FatalErrorIn
        (
            "tmp<tensorField> localAxesRotation::transformTensor"
            "("
                "const tensorField&"
            ")"
        )
            << "tensorField st has different size to tensorField Tr"
            << abort(FatalError);
    }
    return (Rptr_() & tf & Rptr_().T());
}


__HOST____DEVICE__
Foam::tensor Foam::localAxesRotation::transformTensor
(
    const tensor& t
) const
{
    #ifndef __CUDACC__
    notImplemented
    (
        "tensor localAxesRotation::transformTensor(const tensor&) const"
    );
    return tensor::zero;
    #else
    return tensor(0,0,0, 0,0,0, 0,0,0);
    #endif
}


Foam::tmp<Foam::tensorField> Foam::localAxesRotation::transformTensor
(
    const tensorField& tf,
    const labelList& cellMap
) const
{
    if (cellMap.size() != tf.size())
    {
        FatalErrorIn
        (
            "tmp<tensorField> localAxesRotation::transformTensor"
            "("
                "const tensorField&, "
                "const labelList&"
            ")"
        )
            << "tensorField tf has different size to tensorField Tr"
            << abort(FatalError);
    }

    const tensorField& R = Rptr_();
    const tensorField Rtr(R.T());
    tmp<tensorField> tt(new tensorField(cellMap.size()));
    tensorField& t = tt();
    forAll(cellMap, i)
    {
        const label cellI = cellMap[i];
        t[i] = R[cellI] & tf[i] & Rtr[cellI];
    }

    return tt;
}

Foam::tmp<Foam::tensorgpuField> Foam::localAxesRotation::transformTensor
(
    const tensorgpuField& tf
) const
{
    if (Rptr_->size() != tf.size())
    {
        FatalErrorIn
        (
            "tmp<tensorgpuField> localAxesRotation::transformTensor"
            "("
                "const tensorgpuField&"
            ")"
        )
            << "tensorgpuField st has different size to tensorField Tr"
            << abort(FatalError);
    }
    return (Rgpuptr_() & tf & Rgpuptr_().T());
}

namespace Foam
{
	struct localAxesRotationTransform{
		__HOST____DEVICE__
		tensor operator ()(const tensor& t, const thrust::tuple<tensor,tensor>& tr){
			return thrust::get<0>(tr) & t & thrust::get<1>(tr);
		}
	};
}

Foam::tmp<Foam::tensorgpuField> Foam::localAxesRotation::transformTensor
(
    const tensorgpuField& tf,
    const labelgpuList& cellMap
) const
{
    if (cellMap.size() != tf.size())
    {
        FatalErrorIn
        (
            "tmp<tensorField> localAxesRotation::transformTensor"
            "("
                "const tensorgpuField&, "
                "const labelgpuList&"
            ")"
        )
            << "tensorgpuField tf has different size to tensorField Tr"
            << abort(FatalError);
    }

    const tensorgpuField& R = Rgpuptr_();
    const tensorgpuField Rtr(R.T());
    tmp<tensorgpuField> tt(new tensorgpuField(cellMap.size()));
    tensorgpuField& t = tt();
/*
    forAll(cellMap, i)
    {
        const label cellI = cellMap[i];
        t[i] = R[cellI] & tf[i] & Rtr[cellI];
    }
*/
    thrust::transform(tf.begin(),tf.end(),
                      thrust::make_zip_iterator(thrust::make_tuple(
                                                                   thrust::make_permutation_iterator(R.begin(),cellMap.begin()),
                                                                   thrust::make_permutation_iterator(Rtr.begin(),cellMap.begin())
                                                                  )),
                      t.begin(),
                      localAxesRotationTransform());

    return tt;
}

Foam::tmp<Foam::symmTensorField> Foam::localAxesRotation::transformVector
(
    const vectorField& vf
) const
{
    if (Rptr_->size() != vf.size())
    {
        FatalErrorIn("localAxesRotation::transformVector(const vectorField&)")
            << "tensorField vf has different size to tensorField Tr"
            << abort(FatalError);
    }

    tmp<symmTensorField> tfld(new symmTensorField(Rptr_->size()));
    symmTensorField& fld = tfld();

    const tensorField& R = Rptr_();
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R[i], vf[i]);
    }
    return tfld;
}


__HOST____DEVICE__
Foam::symmTensor Foam::localAxesRotation::transformVector
(
    const vector& v
) const
{
    #ifndef __CUDACC__
    notImplemented
    (
        "tensor localAxesRotation::transformVector(const vector&) const"
    );
    return symmTensor::zero;
    #else
    return symmTensor(0,0,0, 0,0,0);
    #endif
}


void Foam::localAxesRotation::write(Ostream& os) const
{
     os.writeKeyword("e3") << e3() << token::END_STATEMENT << nl;
}


// ************************************************************************* //
