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

\*---------------------------------------------------------------------------*/

#include "turbulentInletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct turbulentInletSetupFunctor
{
    __device__
    stateType operator()(const int& id)
    {
        return stateType(id);
    }
};

template<class Type>
struct turbulentInletRandomiseFunctor
{
    __device__
    inline Type operator()(stateType& state);
};

template<>
__device__
scalar turbulentInletRandomiseFunctor<scalar>::operator()(stateType& state)
{
    distributionType dist(0,1);
    return dist(state);
}

template<>
__device__
vector turbulentInletRandomiseFunctor<vector>::operator()(stateType& state)
{
    stateType localState = state;
    distributionType dist(0,1);

    vector out
    (
        dist(localState),
        dist(localState),
        dist(localState)
    );

    state = localState;
    return out;
}

template<>
__device__
tensor turbulentInletRandomiseFunctor<tensor>::operator()(stateType& state)
{
    stateType localState = state;
    distributionType dist(0,1);

    tensor out
    (
        dist(localState),
        dist(localState),
        dist(localState),

        dist(localState),
        dist(localState),
        dist(localState),

        dist(localState),
        dist(localState),
        dist(localState)
    );

    state = localState;
    return out;
}

template<>
__device__
sphericalTensor turbulentInletRandomiseFunctor<sphericalTensor>::operator()(stateType& state)
{
    distributionType dist(0,1);
    return sphericalTensor(dist(state));
}

template<>
__device__
symmTensor turbulentInletRandomiseFunctor<symmTensor>::operator()(stateType& state)
{
    stateType localState = state;
    distributionType dist(0,1);
    symmTensor out
    (
        dist(localState),
        dist(localState),
        dist(localState),

        dist(localState),
        dist(localState),
        dist(localState)
    );

    state = localState;
    return out;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
turbulentInletFvPatchField<Type>::turbulentInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    states_(p.size()),
    fluctuationScale_(pTraits<Type>::zero),
    referenceField_(p.size()),
    alpha_(0.1),
    curTimeIndex_(-1)
{
    updateRandomState();
}


template<class Type>
turbulentInletFvPatchField<Type>::turbulentInletFvPatchField
(
    const turbulentInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    states_(p.size()),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_, mapper),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{
    updateRandomState();
}


template<class Type>
turbulentInletFvPatchField<Type>::turbulentInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    states_(p.size()),
    fluctuationScale_(pTraits<Type>(dict.lookup("fluctuationScale"))),
    referenceField_("referenceField", dict, p.size()),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            gpuField<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(referenceField_);
    }

    updateRandomState();
}


template<class Type>
turbulentInletFvPatchField<Type>::turbulentInletFvPatchField
(
    const turbulentInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    states_(ptf.states_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


template<class Type>
turbulentInletFvPatchField<Type>::turbulentInletFvPatchField
(
    const turbulentInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    states_(ptf.states_),
    fluctuationScale_(ptf.fluctuationScale_),
    referenceField_(ptf.referenceField_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void turbulentInletFvPatchField<Type>::updateRandomState()
{
    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+this->size(),
        states_.begin(),
        turbulentInletSetupFunctor()
    );
}

template<class Type>
void turbulentInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    referenceField_.autoMap(m);
}


template<class Type>
void turbulentInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelgpuList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const turbulentInletFvPatchField<Type>& tiptf =
        refCast<const turbulentInletFvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void turbulentInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        gpuField<Type>& patchField = *this;

        gpuField<Type> randomField(this->size());

        thrust::transform
        (
            states_.begin(),
            states_.end(),
            randomField.begin(),
            turbulentInletRandomiseFunctor<Type>()
        );

        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
        scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

        patchField =
            (1 - alpha_)*patchField
          + alpha_*
            (
                referenceField_
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<Type>::one,
                    fluctuationScale_
                )*mag(referenceField_)
            );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void turbulentInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("fluctuationScale")
        << fluctuationScale_ << token::END_STATEMENT << nl;
    referenceField_.writeEntry("referenceField", os);
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
