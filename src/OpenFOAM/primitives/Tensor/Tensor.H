/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::Tensor

Description
    Templated 3D tensor derived from VectorSpace adding construction from
    9 components, element access using xx(), xy() etc. member functions and
    the inner-product (dot-product) and outer-product of two Vectors
    (tensor-product) operators.

SourceFiles
    TensorI.H

\*---------------------------------------------------------------------------*/

#ifndef Tensor_H
#define Tensor_H

#include "Vector.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Cmpt>
class SymmTensor;

/*---------------------------------------------------------------------------*\
                           Class Tensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor
:
    public VectorSpace<Tensor<Cmpt>, Cmpt, 9>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Tensor<label> labelType;


    // Member constants

        enum
        {
            rank = 2 // Rank of Tensor is 2
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const Tensor zero;
        static const Tensor one;
        static const Tensor max;
        static const Tensor min;
        static const Tensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };


    // Constructors

        //- Construct null
        __HOST____DEVICE__
        inline Tensor();

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        __HOST____DEVICE__
        inline Tensor(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>&);

        //- Construct given SphericalTensor
        __HOST____DEVICE__
        inline Tensor(const SphericalTensor<Cmpt>&);

        //- Construct given SymmTensor
        __HOST____DEVICE__
        inline Tensor(const SymmTensor<Cmpt>&);

        //- Construct given triad
        __HOST____DEVICE__
        inline Tensor(const Vector<Vector<Cmpt> >&);

        //- Construct given the three vector components
        __HOST____DEVICE__
        inline Tensor
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z
        );

        //- Construct given the nine components
        __HOST____DEVICE__
        inline Tensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
            const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
            const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
        );

        //- Construct from Istream
        Tensor(Istream&);


    // Member Functions

        // Access

            __HOST____DEVICE__
            inline const Cmpt& xx() const;
            __HOST____DEVICE__
            inline const Cmpt& xy() const;
            __HOST____DEVICE__
            inline const Cmpt& xz() const;
            __HOST____DEVICE__
            inline const Cmpt& yx() const;
            __HOST____DEVICE__
            inline const Cmpt& yy() const;
            __HOST____DEVICE__
            inline const Cmpt& yz() const;
            __HOST____DEVICE__
            inline const Cmpt& zx() const;
            __HOST____DEVICE__
            inline const Cmpt& zy() const;
            __HOST____DEVICE__
            inline const Cmpt& zz() const;

            __HOST____DEVICE__
            inline Cmpt& xx();
            __HOST____DEVICE__
            inline Cmpt& xy();
            __HOST____DEVICE__
            inline Cmpt& xz();
            __HOST____DEVICE__
            inline Cmpt& yx();
            __HOST____DEVICE__
            inline Cmpt& yy();
            __HOST____DEVICE__
            inline Cmpt& yz();
            __HOST____DEVICE__
            inline Cmpt& zx();
            __HOST____DEVICE__
            inline Cmpt& zy();
            __HOST____DEVICE__
            inline Cmpt& zz();

            // Access vector components.

            __HOST____DEVICE__
            inline Vector<Cmpt> x() const;
            __HOST____DEVICE__
            inline Vector<Cmpt> y() const;
            __HOST____DEVICE__
            inline Vector<Cmpt> z() const;
            __HOST____DEVICE__
            inline Vector<Cmpt> vectorComponent(const direction) const;

        //- Transpose
            __HOST____DEVICE__
        inline Tensor<Cmpt> T() const;


    // Member Operators

        //- Assign to a SphericalTensor
        __HOST____DEVICE__
        inline void operator=(const SphericalTensor<Cmpt>&);

        //- Assign to a SymmTensor
        __HOST____DEVICE__
        inline void operator=(const SymmTensor<Cmpt>&);

        //- Assign to a triad
        __HOST____DEVICE__
        inline void operator=(const Vector<Vector<Cmpt> >&);
};


template<class Cmpt>
class typeOfRank<Cmpt, 2>
{
public:

    typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "TensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
