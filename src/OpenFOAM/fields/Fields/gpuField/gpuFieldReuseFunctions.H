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

#ifndef gpuFieldReuseFunctions_H
#define gpuFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR, class Type1>
class reusegpuTmp
{
public:

    static tmp<gpuField<TypeR> > New(const tmp<gpuField<Type1> >& tf1)
    {
        return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
    }

    static void clear(const tmp<gpuField<Type1> >& tf1)
    {
        tf1.clear();
    }
};


template<class TypeR>
class reusegpuTmp<TypeR, TypeR>
{
public:

    static tmp<gpuField<TypeR> > New(const tmp<gpuField<TypeR> >& tf1)
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
        }
    }

    static void clear(const tmp<gpuField<TypeR> >& tf1)
    {
        if (tf1.isTmp())
        {
            tf1.ptr();
        }
    }
};


template<class TypeR, class Type1, class Type12, class Type2>
class reusegpuTmpTmp
{
public:

    static tmp<gpuField<TypeR> > New
    (
        const tmp<gpuField<Type1> >& tf1,
        const tmp<gpuField<Type2> >& tf2
    )
    {
        return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
    }

    static void clear
    (
        const tmp<gpuField<Type1> >& tf1,
        const tmp<gpuField<Type2> >& tf2
    )
    {
        tf1.clear();
        tf2.clear();
    }
};


template<class TypeR, class Type1, class Type12>
class reusegpuTmpTmp<TypeR, Type1, Type12, TypeR>
{
public:

    static tmp<gpuField<TypeR> > New
    (
        const tmp<gpuField<Type1> >& tf1,
        const tmp<gpuField<TypeR> >& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
        }
    }

    static void clear
    (
        const tmp<gpuField<Type1> >& tf1,
        const tmp<gpuField<TypeR> >& tf2
    )
    {
        tf1.clear();
        if (tf2.isTmp())
        {
            tf2.ptr();
        }
    }
};


template<class TypeR, class Type2>
class reusegpuTmpTmp<TypeR, TypeR, TypeR, Type2>
{
public:

    static tmp<gpuField<TypeR> > New
    (
        const tmp<gpuField<TypeR> >& tf1,
        const tmp<gpuField<Type2> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else
        {
            return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
        }
    }

    static void clear
    (
        const tmp<gpuField<TypeR> >& tf1,
        const tmp<gpuField<Type2> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            tf1.ptr();
        }
        tf2.clear();
    }
};


template<class TypeR>
class reusegpuTmpTmp<TypeR, TypeR, TypeR, TypeR>
{
public:

    static tmp<gpuField<TypeR> > New
    (
        const tmp<gpuField<TypeR> >& tf1,
        const tmp<gpuField<TypeR> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else if (tf2.isTmp())
        {
            return tf2;
        }
        else
        {
            return tmp<gpuField<TypeR> >(new gpuField<TypeR>(tf1().size()));
        }
    }

    static void clear
    (
        const tmp<gpuField<TypeR> >& tf1,
        const tmp<gpuField<TypeR> >& tf2
    )
    {
        if (tf1.isTmp())
        {
            tf1.ptr();
            tf2.clear();
        }
        else if (tf2.isTmp())
        {
            tf1.clear();
            tf2.ptr();
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
