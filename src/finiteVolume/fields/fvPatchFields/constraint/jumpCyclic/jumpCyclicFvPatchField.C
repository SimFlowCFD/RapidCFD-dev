/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "jumpCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const jumpCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict)
{
    // Call this evaluation in derived classes
    //this->evaluate(Pstream::blocking);
}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const jumpCyclicFvPatchField<Type>& ptf
)
:
    cyclicFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const jumpCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
	template<class Type>
	struct jumpCyclicTransformFunctor{
		const tensor t;
		jumpCyclicTransformFunctor(tensor _t): t(_t){}
		__HOST____DEVICE__
		Type operator()(const Type& iF, const Type& jF){
			return transform(t,iF) - jF;
		}
	};
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::jumpCyclicFvPatchField<Type>::patchNeighbourField() const
{
    const gpuField<Type>& iField = this->internalField();
    const labelgpuList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    tmp<gpuField<Type> > tpnf(new gpuField<Type>(this->size()));
    gpuField<Type>& pnf = tpnf();

    gpuField<Type> jf(this->jump());
    if (!this->cyclicPatch().owner())
    {
        jf *= -1.0;
    }

    if (this->doTransform())
    {
        tensor t = this->forwardT().first();
        
        thrust::transform(thrust::make_permutation_iterator(iField.begin(),nbrFaceCells.begin()),
                          thrust::make_permutation_iterator(iField.begin(),nbrFaceCells.begin()+this->size()),
                          jf.begin(),
                          pnf.begin(),
                          jumpCyclicTransformFunctor<Type>(t));
        /*
        forAll(*this, facei)
        {
            pnf[facei] = transform
            (
                this->forwardT()[0], iField[nbrFaceCells[facei]]
            ) - jf[facei];
        }
        */
    }
    else
    {
		thrust::transform(thrust::make_permutation_iterator(iField.begin(),nbrFaceCells.begin()),
		                  thrust::make_permutation_iterator(iField.begin(),nbrFaceCells.end()),
		                  jf.begin(),
		                  pnf.begin(),
		                  subtractOperatorFunctor<Type,Type,Type>());
		/*
        forAll(*this, facei)
        {
            pnf[facei] = iField[nbrFaceCells[facei]] - jf[facei];
        }
        */
    }

    return tpnf;
}


template<class Type>
void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField& psiInternal,
    const scalargpuField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    notImplemented
    (
        "void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrix"
        "("
            "scalarField&, "
            "const scalarField&, "
            "const scalarField&, "
            "const direction, "
            "const Pstream::commsTypes"
        ") const"
    );
}

namespace Foam
{
	template<class Type>
	struct jumpCyclicUpdateInterfaceMatrixFunctor : public std::binary_function<label,Type,Type>{
        const scalar* coeff;
        const Type* pnf;
        const label* neiStart;
        const label* losort;
        jumpCyclicUpdateInterfaceMatrixFunctor(const scalar* _coeff,
                                 const Type* _pnf,
                                 const label* _neiStart,
                                 const label* _losort):
             coeff(_coeff),
             pnf(_pnf),
             neiStart(_neiStart),
             losort(_losort)
        {}
        __HOST____DEVICE__
        Type operator()(const label& id,const Type& s){
            Type out = s;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out -= coeff[face]*pnf[face];
            }

            return out;
        }
    };
}
template<class Type>
void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    gpuField<Type>& result,
    const gpuField<Type>& psiInternal,
    const scalargpuField& coeffs,
    const Pstream::commsTypes
) const
{
    gpuField<Type> pnf(this->size());

    const labelgpuList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    // only apply jump to original field
    if (&psiInternal == &this->internalField())
    {
        gpuField<Type> jf(this->jump());

        if (!this->cyclicPatch().owner())
        {
            jf *= -1.0;
        }
        
        thrust::transform(thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()),
		                  thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()+this->size()),
		                  jf.begin(),
		                  pnf.begin(),
		                  subtractOperatorFunctor<Type,Type,Type>());
        /*
        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]] - jf[facei];
        }
        */
    }
    else
    {
		thrust::copy(thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()),
		             thrust::make_permutation_iterator(psiInternal.begin(),nbrFaceCells.begin()+this->size()),
		             pnf.begin());
		/*
        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]];
        }
        */
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf);

/*
    // Multiply the field by coefficients and add into the result
    const labelgpuList& faceCells = this->cyclicPatch().faceCells();
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
    */
    label patchI = this->patch().index();
    
    const labelgpuList& pcells = this->patch().boundaryMesh().mesh().lduAddr().patchSortCells(patchI);

	const labelgpuList& losort = this->patch().boundaryMesh().mesh().lduAddr().patchSortAddr(patchI);
	const labelgpuList& losortStart = this->patch().boundaryMesh().mesh().lduAddr().patchSortStartAddr(patchI);
    
    thrust::transform(thrust::make_counting_iterator(0),thrust::make_counting_iterator(0)+pcells.size(),
					  thrust::make_permutation_iterator(result.begin(),pcells.begin()),
					  thrust::make_permutation_iterator(result.begin(),pcells.begin()),
                      jumpCyclicUpdateInterfaceMatrixFunctor<Type>(coeffs.data(),
                                               pnf.data(),
                                               losortStart.data(),
                                               losort.data()
                                               ));
}


// ************************************************************************* //
