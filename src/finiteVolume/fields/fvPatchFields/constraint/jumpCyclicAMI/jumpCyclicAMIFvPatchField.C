/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "jumpCyclicAMIFvPatchField.H"
#include "transformField.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::jumpCyclicAMIFvPatchField<Type>::jumpCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMIFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::jumpCyclicAMIFvPatchField<Type>::jumpCyclicAMIFvPatchField
(
    const jumpCyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMIFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::jumpCyclicAMIFvPatchField<Type>::jumpCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMIFvPatchField<Type>(p, iF, dict)
{
    // Call this evaluation in derived classes
    //this->evaluate(Pstream::blocking);
}


template<class Type>
Foam::jumpCyclicAMIFvPatchField<Type>::jumpCyclicAMIFvPatchField
(
    const jumpCyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMIFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::jumpCyclicAMIFvPatchField<Type>::jumpCyclicAMIFvPatchField
(
    const jumpCyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMIFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::jumpCyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const gpuField<Type>& iField = this->internalField();
    const labelgpuList& nbrFaceCells =
        this->cyclicAMIPatch().cyclicAMIPatch().neighbPatch().getFaceCells();

    gpuField<Type> pnf(iField, nbrFaceCells);
    tmp<gpuField<Type> > tpnf;

    if (this->cyclicAMIPatch().applyLowWeightCorrection())
    {
        tpnf =
            this->cyclicAMIPatch().interpolate
            (
                pnf,
                this->patchInternalField()()
            );
    }
    else
    {
        tpnf = this->cyclicAMIPatch().interpolate(pnf);
    }

    if (this->doTransform())
    {
        tpnf = transform(this->getForwardT(), tpnf);
    }

    tmp<gpuField<Type> > tjf = jump();
    if (!this->cyclicAMIPatch().owner())
    {
        tjf = -tjf;
    }

    return tpnf - tjf;
}


template<class Type>
void Foam::jumpCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
        "void Foam::jumpCyclicAMIFvPatchField<Type>::updateInterfaceMatrix"
        "("
            "scalarField&, "
            "const scalarField&, "
            "const scalarField& coeffs,"
            "const direction, "
            "const Pstream::commsTypes"
        ") const"
    );
}


template<class Type>
void Foam::jumpCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    gpuField<Type>& result,
    const gpuField<Type>& psiInternal,
    const scalargpuField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelgpuList& nbrFaceCells =
        this->cyclicAMIPatch().cyclicAMIPatch().neighbPatch().getFaceCells();

    gpuField<Type> pnf(psiInternal, nbrFaceCells);

    if (this->cyclicAMIPatch().applyLowWeightCorrection())
    {
        pnf =
            this->cyclicAMIPatch().interpolate
            (
                pnf,
                this->patchInternalField()()
            );

    }
    else
    {
        pnf = this->cyclicAMIPatch().interpolate(pnf);
    }

    // only apply jump to original field
    if (&psiInternal == &this->internalField())
    {
        gpuField<Type> jf(this->jump());
        if (!this->cyclicAMIPatch().owner())
        {
            jf *= -1.0;
        }

        pnf -= jf;
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf);

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


// ************************************************************************* //
