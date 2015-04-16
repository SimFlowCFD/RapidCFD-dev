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

#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField"
            "("
                "const cyclicAMIFvPatchField<Type>& ,"
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
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
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

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const gpuField<Type>& iField = this->internalField();
    const labelgpuList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().getFaceCells();

    gpuField<Type> pnf(iField, nbrFaceCells);

    tmp<gpuField<Type> > tpnf;
    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf, this->patchInternalField()());
    }
    else
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf);
    }

    if (doTransform())
    {
        tpnf() = transform(getForwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->dimensionedInternalField()
        );

    return refCast<const cyclicAMIFvPatchField<Type> >
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}

template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalargpuField& result,
    const scalargpuField& psiInternal,
    const scalargpuField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelgpuList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().getFaceCells();

    scalargpuField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        scalargpuField pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

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
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    gpuField<Type>& result,
    const gpuField<Type>& psiInternal,
    const scalargpuField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelgpuList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().getFaceCells();

    gpuField<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        gpuField<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

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
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
