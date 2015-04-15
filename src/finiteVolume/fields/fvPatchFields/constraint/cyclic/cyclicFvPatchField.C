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

#include "cyclicFvPatchField.H"
#include "transformField.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isA<cyclicFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField"
            "("
                "const cyclicFvPatchField<Type>& ,"
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isA<cyclicFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField"
            "("
                "const fvPatch&, "
                "const Field<Type>&, "
                "const dictionary&"
            ")",
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    this->evaluate(Pstream::blocking);
}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<gpuField<Type> > cyclicFvPatchField<Type>::patchNeighbourField() const
{
    const gpuField<Type>& iField = this->internalField();
    const labelgpuList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().getFaceCells();

    tmp<gpuField<Type> > tpnf(new gpuField<Type>(this->size()));
    gpuField<Type>& pnf = tpnf();


    if (doTransform())
    {
        tensor t = forwardT()[0];
        
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                iField.begin(),
                nbrFaceCells.begin()
            ),
            thrust::make_permutation_iterator
            (
                iField.begin(),
                nbrFaceCells.end()
            ),
            pnf.begin(),
            transformBinaryFunctionSFFunctor<tensor,Type,Type>(t)
        );
    }
    else
    {
        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                iField.begin(),
                nbrFaceCells.begin()
            ),
            thrust::make_permutation_iterator
            (
                iField.begin(),
                nbrFaceCells.end()
            ),
            pnf.begin()
        );
    }

    return tpnf;
}


template<class Type>
const cyclicFvPatchField<Type>& cyclicFvPatchField<Type>::neighbourPatchField()
const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
    static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
    (
        this->dimensionedInternalField()
    );

    return refCast<const cyclicFvPatchField<Type> >
    (
        fld.boundaryField()[this->cyclicPatch().neighbPatchID()]
    );
}

template<class Type>
struct updateCyclicInterfaceMatrixFunctor
{
   __HOST____DEVICE__
   Type operator()(const Type& s, const thrust::tuple<scalar,Type>& c)
   {
       return s - thrust::get<0>(c) * thrust::get<1>(c);
   }
};

template<class Type>
void cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField& psiInternal,
    const scalargpuField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelgpuList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().getFaceCells();

    scalargpuField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    matrixPatchOperation
    (
        this->patch().index(),
        result,
        this->patch().boundaryMesh().mesh().lduAddr(),
        matrixInterfaceFunctor<scalar>
        (
            coeffs.data(),
            pnf.data()
        )
    );
}


template<class Type>
void cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    gpuField<Type>& result,
    const gpuField<Type>& psiInternal,
    const scalargpuField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelgpuList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().getFaceCells();

    gpuField<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    // Multiply the field by coefficients and add into the result
    matrixPatchOperation
    (
        this->patch().index(),
        result,
        this->patch().boundaryMesh().mesh().lduAddr(),
        matrixInterfaceFunctor<Type>
        (
            coeffs.data(),
            pnf.data()
        )
    );
}


template<class Type>
void cyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
