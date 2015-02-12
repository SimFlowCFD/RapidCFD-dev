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
    Foam::gpuFieldMapper

Description
    Abstract base class to hold the gpuField mapping addressing and weights.

\*---------------------------------------------------------------------------*/

#ifndef gpuFieldMapper_H
#define gpuFieldMapper_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gpuFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class gpuFieldMapper
{

public:

    // Constructors

        //- Null constructor
        gpuFieldMapper()
        {}


    //- Destructor
    virtual ~gpuFieldMapper()
    {}


    // Member Functions

        virtual label size() const = 0;

        virtual bool direct() const = 0;

        //- Are there unmapped values? I.e. do all size() elements get
        //  get value
        virtual bool hasUnmapped() const = 0;

        virtual const labelgpuList& directAddressing() const
        {
            FatalErrorIn("gpuFieldMapper::directAddressing() const")
                << "attempt to access null direct addressing"
                << abort(FatalError);

            return labelgpuList::null();
        }

        virtual const labelgpuListList& addressing() const
        {
            FatalErrorIn("gpuFieldMapper::addressing() const")
                << "attempt to access null interpolation addressing"
                << abort(FatalError);

            return labelgpuListList::null();
        }

        virtual const scalargpuListList& weights() const
        {
            FatalErrorIn("gpuFieldMapper::weights() const")
                << "attempt to access null interpolation weights"
                << abort(FatalError);

            return scalargpuListList::null();
        }


    // Member Operators

        template<class Type>
        tmp<gpuField<Type> > operator()(const gpuField<Type>& f) const
        {
            return tmp<gpuField<Type> >(new gpuField<Type>(f, *this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
