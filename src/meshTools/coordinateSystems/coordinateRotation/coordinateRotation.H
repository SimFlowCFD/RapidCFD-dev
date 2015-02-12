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

Class
    Foam::coordinateRotation

Description
    Abstract base class for coordinate rotation

    \verbatim
        coordinateRotation
        {
            type        axesRotation
            e1          (1 0 0);
            e2          (0 1 0);
        }
    \endverbatim

    Types of coordinateRotation:
    1) axesRotation
    2) STARCDRotation
    3) localAxesRotation
    4) EulerCoordinateRotation



\*---------------------------------------------------------------------------*/

#ifndef coordinateRotation_H
#define coordinateRotation_H

#include "vector.H"
#include "tensor.H"
#include "tensorField.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "objectRegistry.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class coordinateRotation
{

public:


    //- Runtime type information
    TypeName("coordinateRotation");


    // Declare run-time constructor selection table
    // for constructors with dictionary and objectRegistry
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateRotation,
        objectRegistry,
        (
            const dictionary& dict, const objectRegistry& obr
        ),
        (dict, obr)
    );


    // Declare run-time constructor selection table
    // for constructors with dictionary
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateRotation,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

        //- Transform principal
        __HOST____DEVICE__
        symmTensor transformPrincipal(const tensor&, const vector&) const;

    // Selectors

        //- Select constructed from dictionary and objectRegistry
        static autoPtr<coordinateRotation> New
        (
            const dictionary& dict, const objectRegistry& obr
        );

        //- Select constructed from dictionary
        static autoPtr<coordinateRotation> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~coordinateRotation()
    {}


    // Member Functions

        //- Reset rotation to an identity rotation
        virtual void clear() = 0;

        //- Update the rotation for a list of cells
        virtual void updateCells
        (
            const polyMesh& mesh,
            const labelgpuList& cells
        ) = 0;

        //- Return local-to-global transformation tensor
        __HOST____DEVICE__
        virtual const tensor& R() const = 0;

        //- Return global-to-local transformation tensor
        __HOST____DEVICE__
        virtual const tensor& Rtr() const = 0;

        //- Return local Cartesian x-axis
        __HOST____DEVICE__
        virtual const vector e1() const = 0;

        //- Return local Cartesian y-axis
        __HOST____DEVICE__
        virtual const vector e2() const = 0;

        //- Return local Cartesian z-axis
        __HOST____DEVICE__
        virtual const vector e3() const = 0;

        //- Return local-to-global transformation tensor
        virtual const tensorField& Tr() const = 0;

        //- Return true if the rotation tensor is uniform
        __HOST____DEVICE__
        virtual bool uniform() const
        {
            return true;
        }

        //- Transform vectorField using transformation tensor field
        virtual tmp<vectorField> transform(const vectorField& st) const = 0;

        //- Transform vectorgpuField using transformation tensor field
        virtual tmp<vectorgpuField> transform(const vectorgpuField& st) const = 0;

        //- Transform vector using transformation tensor
        __HOST____DEVICE__
        virtual vector transform(const vector& st) const = 0;

        //- Inverse transform vectorField using transformation tensor field
        virtual tmp<vectorField> invTransform(const vectorField& st) const = 0;

        //- Inverse transform vectorgpuField using transformation tensor field
        virtual tmp<vectorgpuField> invTransform(const vectorgpuField& st) const = 0;

        //- Inverse transform vector using transformation tensor
        __HOST____DEVICE__
        virtual vector invTransform(const vector& st) const = 0;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorField> transformTensor
        (
            const tensorField& st
        ) const = 0;

        //- Transform tensor sub-field using transformation tensorField
        virtual tmp<tensorField> transformTensor
        (
            const tensorField& st,
            const labelList& cellMap
        ) const = 0;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorgpuField> transformTensor
        (
            const tensorgpuField& st
        ) const = 0;

        //- Transform tensor sub-field using transformation tensorField
        virtual tmp<tensorgpuField> transformTensor
        (
            const tensorgpuField& st,
            const labelgpuList& cellMap
        ) const = 0;

        //- Transform tensor using transformation tensorField
        __HOST____DEVICE__
        virtual tensor transformTensor(const tensor& st) const = 0;

        //- Transform vectorField using transformation tensorField and return
        // symmetrical tensorField
        virtual tmp<symmTensorField> transformVector
        (
            const vectorField& st
        ) const = 0;

        virtual tmp<symmTensorgpuField> transformVector
        (
            const vectorgpuField& st
        ) const = 0;

        //- Transform vector using transformation tensor and return
        // symmetrical tensor
        __HOST____DEVICE__
        virtual symmTensor transformVector(const vector& st) const = 0;


    // Write

            //- Write
            virtual void write(Ostream&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "coordinateRotationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
