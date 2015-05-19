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

Class
    Foam::mixedFixedValueSlipFvPatchField

Description
    A mixed boundary type that blends between fixedValue and slip, as opposed
    to the standard mixed condition that blends between fixedValue and
    fixedGradient; required to implement maxwellSlipU condition.

SourceFiles
    mixedFixedValueSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedFixedValueSlipFvPatchField_H
#define mixedFixedValueSlipFvPatchField_H

#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixedFixedValueSlipFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedFixedValueSlipFvPatchField
:
    public transformFvPatchField<Type>
{
    // Private data

        //- Value field used for boundary condition
        gpuField<Type> refValue_;

        //- Fraction (0-1) of value used for boundary condition
        scalargpuField valueFraction_;

public:

    //- Runtime type information
    TypeName("mixedFixedValueSlip");


    // Constructors

        //- Construct from patch and internal field
        mixedFixedValueSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixedFixedValueSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedFixedValueSlipFvPatchField
        //- onto a new patch
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedFixedValueSlipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedFixedValueSlipFvPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelgpuList&
            );

        // Return defining fields

            virtual gpuField<Type>& refValue()
            {
                return refValue_;
            }

            virtual const gpuField<Type>& refValue() const
            {
                return refValue_;
            }

            virtual scalargpuField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalargpuField& valueFraction() const
            {
                return valueFraction_;
            }

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<gpuField<Type> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return face-gradient transform diagonal
            virtual tmp<gpuField<Type> > snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const gpuList<Type>&) {}
        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const gpuField<Type>&) {}
        virtual void operator-=(const gpuField<Type>&) {}

        virtual void operator*=(const gpuField<scalar>&) {}
        virtual void operator/=(const gpuField<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixedFixedValueSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
