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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
inline const DimensionedField<Type, GeoMesh>&
DimensionedField<Type, GeoMesh>::null()
{
    return *reinterpret_cast< DimensionedField<Type, GeoMesh>* >(0);
}


template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
DimensionedField<Type, GeoMesh>::mesh() const
{
    return mesh_;
}


template<class Type, class GeoMesh>
inline const dimensionSet& DimensionedField<Type, GeoMesh>::dimensions() const
{
    return dimensions_;
}

template<class Type, class GeoMesh>
inline dimensionSet& DimensionedField<Type, GeoMesh>::dimensions()
{
    return dimensions_;
}

/*
template<class Type, class GeoMesh>
inline Field<Type>& DimensionedField<Type, GeoMesh>::field()
{
    return *this;
}
*/

template<class Type, class GeoMesh>
inline const gpuField<Type>& DimensionedField<Type, GeoMesh>::getField() const
{
    return field_;
}

template<class Type, class GeoMesh>
inline gpuField<Type>& DimensionedField<Type, GeoMesh>::getField()
{
    return field_;
}

template<class Type, class GeoMesh>
inline DimensionedField<Type, GeoMesh>::operator const gpuField<Type>&() const
{
    return field_;
}

template<class Type, class GeoMesh>
inline DimensionedField<Type, GeoMesh>::operator gpuField<Type>&()
{
    return field_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
