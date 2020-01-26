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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "IOstreams.H"
#include "Switch.H"
#include "BasicCache.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduMatrix, 1);

    class lduMatrixCache
    {
        static PtrList<scalargpuField> diagCache;
        static PtrList<scalargpuField> lowerCache;
        static PtrList<scalargpuField> upperCache;

        static PtrList<scalargpuField> lowerSortCache;
        static PtrList<scalargpuField> upperSortCache;

        public:

        static scalargpuField& diag(label level, label size)
        {
            return cache::retrieve(diagCache,level,size);
        }

        static scalargpuField& lower(label level, label size)
        {
            return cache::retrieve(lowerCache,level,size);
        }

        static scalargpuField& upper(label level, label size)
        {
            return cache::retrieve(upperCache,level,size);
        }

        static scalargpuField* lowerSort(label level, label size)
        {
            return &cache::retrieve(lowerSortCache,level,size);
        }

        static scalargpuField* upperSort(label level, label size)
        {
            return &cache::retrieve(upperSortCache,level,size);
        }
    };

    PtrList<scalargpuField> lduMatrixCache::diagCache(1);
    PtrList<scalargpuField> lduMatrixCache::lowerCache(1);
    PtrList<scalargpuField> lduMatrixCache::upperCache(1);
    PtrList<scalargpuField> lduMatrixCache::lowerSortCache(1);
    PtrList<scalargpuField> lduMatrixCache::upperSortCache(1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::lduMatrix::lduMatrix(const lduMesh& mesh)
:
    lduMesh_(mesh),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerSortPtr_(NULL),
    upperSortPtr_(NULL),
    coarsestLevel_(false)
{}


Foam::lduMatrix::lduMatrix(const lduMatrix& A)
:
    lduMesh_(A.lduMesh_),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerSortPtr_(NULL),
    upperSortPtr_(NULL),
    coarsestLevel_(false)
{
    if (A.lowerPtr_)
    {
        lowerPtr_ = new scalargpuField(*(A.lowerPtr_));
    }

    if (A.diagPtr_)
    {
        diagPtr_ = new scalargpuField(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = new scalargpuField(*(A.upperPtr_));
    }
}


Foam::lduMatrix::lduMatrix(lduMatrix& A, bool reUse)
:
    lduMesh_(A.lduMesh_),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerSortPtr_(NULL),
    upperSortPtr_(NULL),
    coarsestLevel_(false)
{
    if (reUse)
    {
        if (A.lowerPtr_)
        {
            lowerPtr_ = A.lowerPtr_;
            A.lowerPtr_ = NULL;
        }

        if (A.diagPtr_)
        {
            diagPtr_ = A.diagPtr_;
            A.diagPtr_ = NULL;
        }

        if (A.upperPtr_)
        {
            upperPtr_ = A.upperPtr_;
            A.upperPtr_ = NULL;
        }
    }
    else
    {
        if (A.lowerPtr_)
        {
            lowerPtr_ = new scalargpuField(*(A.lowerPtr_));
        }

        if (A.diagPtr_)
        {
            diagPtr_ = new scalargpuField(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = new scalargpuField(*(A.upperPtr_));
        }
    }
}


Foam::lduMatrix::lduMatrix(const lduMesh& mesh, Istream& is)
:
    lduMesh_(mesh),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    lowerSortPtr_(NULL),
    upperSortPtr_(NULL),
    coarsestLevel_(false)
{
    Switch hasLow(is);
    Switch hasDiag(is);
    Switch hasUp(is);

    if (hasLow)
    {
        lowerPtr_ = new scalargpuField(is);
    }
    if (hasDiag)
    {
        diagPtr_ = new scalargpuField(is);
    }
    if (hasUp)
    {
        upperPtr_ = new scalargpuField(is);
    }
}


Foam::lduMatrix::~lduMatrix()
{
    if (lowerPtr_)
    {
        delete lowerPtr_;
    }

    if (diagPtr_)
    {
        delete diagPtr_;
    }

    if (upperPtr_)
    {
        delete upperPtr_;
    }
}


Foam::scalargpuField& Foam::lduMatrix::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = new scalargpuField(*upperPtr_);
        }
        else
        {
            lowerPtr_ = new scalargpuField(lduAddr().lowerAddr().size(), 0.0);
        }
    }

    lowerSortPtr_ = NULL;

    return *lowerPtr_;
}


Foam::scalargpuField& Foam::lduMatrix::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ = new scalargpuField(lduAddr().size(), 0.0);
    }

    return *diagPtr_;
}


Foam::scalargpuField& Foam::lduMatrix::upper()
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = new scalargpuField(*lowerPtr_);
        }
        else
        {
            upperPtr_ = new scalargpuField(lduAddr().lowerAddr().size(), 0.0);
        }
    }

    upperSortPtr_ = NULL;

    return *upperPtr_;
}


Foam::scalargpuField& Foam::lduMatrix::lower(const label nCoeffs)
{
    if (!lowerPtr_)
    {
        lowerPtr_ = new scalargpuField(const_cast<const scalargpuField&>(lduMatrixCache::lower(level(),nCoeffs)),nCoeffs);
        if (upperPtr_)
        {
            *lowerPtr_ = *upperPtr_;
        }
        else
        {
            *lowerPtr_ = 0.0;
        }
    }

    lowerSortPtr_ = NULL;

    return *lowerPtr_;
}


Foam::scalargpuField& Foam::lduMatrix::diag(const label size)
{
    if (!diagPtr_)
    {
        diagPtr_ = new scalargpuField(const_cast<const scalargpuField&>(lduMatrixCache::diag(level(),size)),size);
        *diagPtr_ = 0.0;
    }

    return *diagPtr_;
}


Foam::scalargpuField& Foam::lduMatrix::upper(const label nCoeffs)
{
    if (!upperPtr_)
    {
        upperPtr_ = new scalargpuField(const_cast<const scalargpuField&>(lduMatrixCache::upper(level(),nCoeffs)),nCoeffs);

        if (lowerPtr_)
        {
            *upperPtr_ = *lowerPtr_;
        }
        else
        {
            
            *upperPtr_ = 0.0;
        }
    }

    upperSortPtr_ = NULL;

    return *upperPtr_;
}


const Foam::scalargpuField& Foam::lduMatrix::lower() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::lower() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if (lowerPtr_)
    {
        return *lowerPtr_;
    }
    else
    {
        return *upperPtr_;
    }
}


const Foam::scalargpuField& Foam::lduMatrix::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorIn("const scalargpuField& lduMatrix::diag() const")
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


const Foam::scalargpuField& Foam::lduMatrix::upper() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::upper() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if (upperPtr_)
    {
        return *upperPtr_;
    }
    else
    {
        return *lowerPtr_;
    }
}

void Foam::lduMatrix::calcSortCoeffs
(
    scalargpuField& out, 
    const scalargpuField& in
) const
{
    const labelgpuList& sort = lduAddr().losortAddr();

    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            in.begin(),
            sort.begin()
        ),
        thrust::make_permutation_iterator
        (
            in.begin(),
            sort.end()
        ),
        out.begin()
    );
}

const Foam::scalargpuField& Foam::lduMatrix::lowerSort() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::lowerSort() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if(lowerSortPtr_)
    {
        return *lowerSortPtr_;
    }

    if (lowerPtr_)
    {   
        lowerSortPtr_ = lduMatrixCache::lowerSort(level(),lowerPtr_->size());

        calcSortCoeffs(*lowerSortPtr_,*lowerPtr_);
      
        return *lowerSortPtr_;
    }
    else
    {
        if( ! upperSortPtr_)
        {
            upperSortPtr_ = lduMatrixCache::upperSort(level(),upperPtr_->size());

            calcSortCoeffs(*upperSortPtr_,*upperPtr_);
        }
      
        return *upperSortPtr_;
    }
}

const Foam::scalargpuField& Foam::lduMatrix::upperSort() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::upperSort() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if(upperSortPtr_)
    {
        return *upperSortPtr_;
    }

    if (upperPtr_)
    {
        upperSortPtr_ = lduMatrixCache::upperSort(level(),upperPtr_->size());

        calcSortCoeffs(*upperSortPtr_,*upperPtr_);
      
        return *upperSortPtr_;
    }
    else
    {
        if( ! lowerSortPtr_)
        {
            lowerSortPtr_ = lduMatrixCache::lowerSort(level(),lowerPtr_->size());

            calcSortCoeffs(*lowerSortPtr_,*lowerPtr_);
        }
      
        return *lowerSortPtr_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const lduMatrix& ldum)
{
    Switch hasLow = ldum.hasLower();
    Switch hasDiag = ldum.hasDiag();
    Switch hasUp = ldum.hasUpper();

    os  << hasLow << token::SPACE << hasDiag << token::SPACE
        << hasUp << token::SPACE;

    if (hasLow)
    {
        os  << ldum.lower();
    }

    if (hasDiag)
    {
        os  << ldum.diag();
    }

    if (hasUp)
    {
        os  << ldum.upper();
    }

    os.check("Ostream& operator<<(Ostream&, const lduMatrix&");

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<lduMatrix>& ip)
{
    const lduMatrix& ldum = ip.t_;

    Switch hasLow = ldum.hasLower();
    Switch hasDiag = ldum.hasDiag();
    Switch hasUp = ldum.hasUpper();

    os  << "Lower:" << hasLow
        << " Diag:" << hasDiag
        << " Upper:" << hasUp << endl;

    if (hasLow)
    {
        os  << "lower:" << ldum.lower().size() << endl;
    }
    if (hasDiag)
    {
        os  << "diag :" << ldum.diag().size() << endl;
    }
    if (hasUp)
    {
        os  << "upper:" << ldum.upper().size() << endl;
    }


    //if (hasLow)
    //{
    //    os  << "lower contents:" << endl;
    //    forAll(ldum.lower(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.lower()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasDiag)
    //{
    //    os  << "diag contents:" << endl;
    //    forAll(ldum.diag(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.diag()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasUp)
    //{
    //    os  << "upper contents:" << endl;
    //    forAll(ldum.upper(), i)
    //    {
    //        os  << "i:" << i << "\t" << ldum.upper()[i] << endl;
    //    }
    //    os  << endl;
    //}

    os.check("Ostream& operator<<(Ostream&, const lduMatrix&");

    return os;
}


// ************************************************************************* //
