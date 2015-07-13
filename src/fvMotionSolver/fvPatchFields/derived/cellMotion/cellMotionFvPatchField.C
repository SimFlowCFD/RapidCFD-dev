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

#include "cellMotionFvPatchField.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "pointFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF)
{
    fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
    const cellMotionFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
template<class Type>
struct cellMotionFvPatchFieldUpdateCoeffsFunctor
{
    const label* facePoints;
    const point* points;
    const Type* field;
    const Type zero;

    cellMotionFvPatchFieldUpdateCoeffsFunctor
    (
        const label* _facePoints,
        const point* _points,
        const Type* _field
    ):
        facePoints(_facePoints),
        points(_points), 
        field(_field),
        zero(pTraits<Type>::zero)
    {}

    __host__ __device__
    Type operator()(const faceData& face)
    {
        return face.average(facePoints,points,field,zero);
    }
};
}

template<class Type>
void Foam::cellMotionFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& p = this->patch();
    const polyPatch& pp = p.patch();
    const fvMesh& mesh = this->dimensionedInternalField().mesh();
    const pointgpuField& points = mesh.getPoints();

    const faceDatagpuList& faces = pp.getFaces();
    const labelgpuList& facePoints = pp.getFaceNodes();

    word pfName = this->dimensionedInternalField().name();
    pfName.replace("cell", "point");

    const GeometricField<Type, pointPatchField, pointMesh>& pointMotion =
        this->db().objectRegistry::template
            lookupObject<GeometricField<Type, pointPatchField, pointMesh> >
            (pfName);

    thrust::transform
    (
        faces.begin(),
        faces.end(),
        this->begin(),
        cellMotionFvPatchFieldUpdateCoeffsFunctor<Type>
        (
            facePoints.data(),
            points.data(),
            pointMotion.getField().data()
        )
    );

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::cellMotionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}

// ************************************************************************* //
