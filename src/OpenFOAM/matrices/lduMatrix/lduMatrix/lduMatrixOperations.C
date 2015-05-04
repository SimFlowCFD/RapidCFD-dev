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

Description
    lduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "FieldM.H"
#include "lduMatrixSolutionCache.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::lduMatrix::sumDiag()
{
    const scalargpuField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalargpuField& Upper = const_cast<const lduMatrix&>(*this).upper();

    matrixOperation
    (
        diag().begin(),
        diag(),
        lduAddr(),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            Lower.data(),
            unityOp<scalar>()
        ),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            Upper.data(),
            unityOp<scalar>()
        )
    );                    
}

void Foam::lduMatrix::negSumDiag()
{
    const scalargpuField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalargpuField& Upper = const_cast<const lduMatrix&>(*this).upper();

    matrixOperation
    (
        diag().begin(),
        diag(),
        lduAddr(),
        matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >
        (
            Lower.data(),
            negateUnaryOperatorFunctor<scalar,scalar>()
        ),
        matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >
        (
            Upper.data(),
            negateUnaryOperatorFunctor<scalar,scalar>()
        )
    );  
}

void Foam::lduMatrix::sumMagOffDiag
(
    scalargpuField& sumOff
) const
{
    const scalargpuField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalargpuField& Upper = const_cast<const lduMatrix&>(*this).upper();

    matrixOperation
    (
        sumOff.begin(),
        sumOff,
        lduAddr(),
        matrixCoeffsFunctor<scalar,magUnaryFunctionFunctor<scalar,scalar> >
        (
            Upper.data(),
            magUnaryFunctionFunctor<scalar,scalar>()
        ),
        matrixCoeffsFunctor<scalar,magUnaryFunctionFunctor<scalar,scalar> >
        (
            Lower.data(),
            magUnaryFunctionFunctor<scalar,scalar>()
        )
    ); 
}

#define H_FUNCTION_CALL(functionName)                                                       \
functionName                                                                                \
(                                                                                           \
    Hpsi.begin(),                                                                           \
    Hpsi,                                                                                   \
    lduAddr(),                                                                              \
    matrixCoeffsMultiplyFunctor<scalar,scalar,negateUnaryOperatorFunctor<scalar,scalar> >   \
    (                                                                                       \
        psi.data(),                                                                         \
        Upper.data(),                                                                       \
        u.data(),                                                                           \
        negateUnaryOperatorFunctor<scalar,scalar>()                                         \
    ),                                                                                      \
    matrixCoeffsMultiplyFunctor<scalar,scalar,negateUnaryOperatorFunctor<scalar,scalar> >   \
    (                                                                                       \
        psi.data(),                                                                         \
        Lower.data(),                                                                       \
        l.data(),                                                                           \
        negateUnaryOperatorFunctor<scalar,scalar>()                                         \
    )                                                                                       \
);

template<>
void Foam::lduMatrix::H(Foam::gpuField<Foam::scalar>& Hpsi,const Foam::gpuField<Foam::scalar>& psi) const
{
    Hpsi = 0;

    if (lowerPtr_ || upperPtr_)
    {
        bool fastPath = lduMatrixSolutionCache::favourSpeed;

        const scalargpuField& Lower = fastPath?this->lowerSort():this->lower();
        const scalargpuField& Upper = this->upper();

        const labelgpuList& l = fastPath?lduAddr().ownerSortAddr():lduAddr().lowerAddr();
        const labelgpuList& u = lduAddr().upperAddr();
        
        if(fastPath)
        {
            H_FUNCTION_CALL(matrixFastOperation);
        }
        else
        {
            H_FUNCTION_CALL(matrixOperation);
        }                              
    }
}

template<>
Foam::tmp<Foam::gpuField<Foam::scalar> > Foam::lduMatrix::H(const Foam::gpuField<Foam::scalar>& psi) const
{
    tmp<gpuField<scalar> > tHpsi
    (
        new gpuField<scalar>(lduAddr().size(), 0)
    );

    H(tHpsi(),psi);

    return tHpsi;
}

#undef H_FUNCTION_CALL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::operator=(const lduMatrix& A)
{
    if (this == &A)
    {
        FatalError
            << "lduMatrix::operator=(const lduMatrix&) : "
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (A.lowerPtr_)
    {
        lower() = A.lower();
    }
    else if (lowerPtr_)
    {
        delete lowerPtr_;
        lowerPtr_ = NULL;
    }

    if (A.upperPtr_)
    {
        upper() = A.upper();
    }
    else if (upperPtr_)
    {
        delete upperPtr_;
        upperPtr_ = NULL;
    }

    if (A.diagPtr_)
    {
        diag() = A.diag();
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


void Foam::lduMatrix::negate()
{
    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


void Foam::lduMatrix::operator+=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() += A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() += A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningIn("lduMatrix::operator+=(const lduMatrix& A)")
                << "Unknown matrix type combination" << nl
                << "    this :"
                << " diagonal:" << diagonal()
                << " symmetric:" << symmetric()
                << " asymmetric:" << asymmetric() << nl
                << "    A    :"
                << " diagonal:" << A.diagonal()
                << " symmetric:" << A.symmetric()
                << " asymmetric:" << A.asymmetric()
                << endl;
        }
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


void Foam::lduMatrix::operator-=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() -= A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = -A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        if (debug > 1)
        {
            WarningIn("lduMatrix::operator-=(const lduMatrix& A)")
                << "Unknown matrix type combination" << nl
                << "    this :"
                << " diagonal:" << diagonal()
                << " symmetric:" << symmetric()
                << " asymmetric:" << asymmetric() << nl
                << "    A    :"
                << " diagonal:" << A.diagonal()
                << " symmetric:" << A.symmetric()
                << " asymmetric:" << A.asymmetric()
                << endl;
        }
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


void Foam::lduMatrix::operator*=(const scalargpuField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        scalargpuField& upper = *upperPtr_;

        const labelgpuList& l = lduAddr().lowerAddr();

        thrust::transform
        (
            upper.begin(),
            upper.end(),
            thrust::make_permutation_iterator(sf.begin(),l.begin()),
            upper.begin(),
            multiplyOperatorFunctor<scalar,scalar,scalar>()
        );
    }

    if (lowerPtr_)
    {
        scalargpuField& lower = *lowerPtr_;

        const labelgpuList& u = lduAddr().upperAddr();

        thrust::transform
        (
            lower.begin(),
            lower.end(),
            thrust::make_permutation_iterator(sf.begin(),u.begin()),
            lower.begin(),
            multiplyOperatorFunctor<scalar,scalar,scalar>()
        );
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


void Foam::lduMatrix::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }

    upperSortPtr_ = NULL;
    lowerSortPtr_ = NULL;
}


// ************************************************************************* //
