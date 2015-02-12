/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct faceHLduMatrixFunctor
{
    template<class Type,class Tuple>
    __HOST____DEVICE__
    Type operator()(const Tuple& t)
    {
        return thrust::get<0>(t)*thrust::get<1>(t) 
               - thrust::get<2>(t)*thrust::get<3>(t);
    }
};

}

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumDiag()
{
    const gpuField<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const gpuField<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();
    gpuField<DType>& Diag = diag();

    matrixOperation
    (
        Diag.begin(),
        Diag,
        lduAddr(),
        matrixCoeffsFunctor<LUType,unityOp<LUType> >
        (
            Lower.data(),
            unityOp<LUType>()
        ),
        matrixCoeffsFunctor<LUType,unityOp<LUType> >
        (
            Upper.data(),
            unityOp<LUType>()
        )
    );   
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::negSumDiag()
{
    const gpuField<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const gpuField<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();
    gpuField<DType>& Diag = diag();

    matrixOperation
    (
        Diag.begin(),
        Diag,
        lduAddr(),
        matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<LUType,LUType> >
        (
            Lower.data(),
            negateUnaryOperatorFunctor<LUType,LUType>()
        ),
        matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<LUType,LUType> >
        (
            Upper.data(),
            negateUnaryOperatorFunctor<LUType,LUType>()
        )
    );  
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumMagOffDiag
(
    gpuField<LUType>& sumOff
) const
{
    const gpuField<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const gpuField<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();

    matrixOperation
    (
        sumOff.begin(),
        sumOff,
        lduAddr(),
        matrixCoeffsFunctor<scalar,cmptMagUnaryFunctionFunctor<LUType,LUType> >
        (
            Upper.data(),
            cmptMagUnaryFunctionFunctor<LUType,LUType>()
        ),
        matrixCoeffsFunctor<scalar,cmptMagUnaryFunctionFunctor<LUType,LUType> >
        (
            Lower.data(),
            cmptMagUnaryFunctionFunctor<LUType,LUType>()
        )
    ); 
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::gpuField<Type> >
Foam::LduMatrix<Type, DType, LUType>::H(const gpuField<Type>& psi) const
{
    tmp<gpuField<Type> > tHpsi
    (
        new gpuField<Type>(lduAddr().size(), pTraits<Type>::zero)
    );

    if (lowerPtr_ || upperPtr_)
    {
        gpuField<Type> & Hpsi = tHpsi();

        const labelgpuList& l = lduAddr().lowerAddr();
        const labelgpuList& u = lduAddr().upperAddr();

        const gpuField<LUType>& Lower = lower(); 
        const gpuField<LUType>& Upper = upper(); 

        matrixOperation
        (
            Hpsi.begin(),
            Hpsi,
            lduAddr(),
            matrixCoeffsMultiplyFunctor<Type,LUType,negateUnaryOperatorFunctor<Type,Type> >
            (
                psi.data(),
                Upper.data(),
                u.data(),
                negateUnaryOperatorFunctor<Type,Type>()
            ),
            matrixCoeffsMultiplyFunctor<Type,LUType,negateUnaryOperatorFunctor<Type,Type> >
            (
                psi.data(),
                Lower.data(),
                l.data(),
                negateUnaryOperatorFunctor<Type,Type>()
            )
        );    

    }

    return tHpsi;
}

template<class Type, class DType, class LUType>
Foam::tmp<Foam::gpuField<Type> >
Foam::LduMatrix<Type, DType, LUType>::H(const tmp<gpuField<Type> >& tpsi) const
{
    tmp<gpuField<Type> > tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::gpuField<Type> >
Foam::LduMatrix<Type, DType, LUType>::faceH(const gpuField<Type>& psi) const
{
    const gpuField<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const gpuField<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();

    // Take refereces to addressing
    const labelgpuList& l = lduAddr().lowerAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    tmp<gpuField<Type> > tfaceHpsi(new gpuField<Type> (Lower.size()));
    gpuField<Type> & faceHpsi = tfaceHpsi();

    thrust::transform
    (
        thrust::make_zip_iterator(thrust::make_tuple
        (
            Upper.begin(),
            thrust::make_permutation_iterator(psi.begin(),u.begin()),
            Lower.begin(),
            thrust::make_permutation_iterator(psi.begin(),l.begin())
        )),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            Upper.end(),
            thrust::make_permutation_iterator(psi.begin(),u.end()),
            Lower.end(),
            thrust::make_permutation_iterator(psi.begin(),l.end())
        )),
        faceHpsi.begin(),
        faceHLduMatrixFunctor()
    );

    return tfaceHpsi;
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::gpuField<Type> >
Foam::LduMatrix<Type, DType, LUType>::faceH(const tmp<gpuField<Type> >& tpsi) const
{
    tmp<gpuField<Type> > tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator=(const LduMatrix& A)
{
    if (this == &A)
    {
        FatalErrorIn
        (
            "LduMatrix<Type, DType, LUType>::operator=(const LduMatrix&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    if (A.diagPtr_)
    {
        diag() = A.diag();
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

    if (A.lowerPtr_)
    {
        lower() = A.lower();
    }
    else if (lowerPtr_)
    {
        delete lowerPtr_;
        lowerPtr_ = NULL;
    }

    if (A.sourcePtr_)
    {
        source() = A.source();
    }

    interfacesUpper_ = A.interfacesUpper_;
    interfacesLower_ = A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::negate()
{
    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (sourcePtr_)
    {
        sourcePtr_->negate();
    }

    Foam::negate(interfacesUpper_);
    Foam::negate(interfacesLower_);
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator+=(const LduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() += A.diag();
    }

    if (A.sourcePtr_)
    {
        source() += A.source();
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
        FatalErrorIn
        (
            "LduMatrix<Type, DType, LUType>::operator+=(const LduMatrix& A)"
        )   << "Unknown matrix type combination"
            << abort(FatalError);
    }

    interfacesUpper_ += A.interfacesUpper_;
    interfacesLower_ += A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator-=(const LduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    if (A.sourcePtr_)
    {
        source() -= A.source();
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
        FatalErrorIn
        (
            "LduMatrix<Type, DType, LUType>::operator-=(const LduMatrix& A)"
        )   << "Unknown matrix type combination"
            << abort(FatalError);
    }

    interfacesUpper_ -= A.interfacesUpper_;
    interfacesLower_ -= A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator*=
(
    const scalargpuField& sf
)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (sourcePtr_)
    {
        *sourcePtr_ *= sf;
    }

    if (upperPtr_)
    {
        gpuField<LUType>& upper = *upperPtr_;

        const labelgpuList& l = lduAddr().lowerAddr();

        thrust::transform
        (
            upper.begin(),
            upper.end(),
            thrust::make_permutation_iterator(sf.begin(),l.begin()),
            upper.begin(),
            multiplyOperatorFunctor<LUType,LUType,LUType>()
        );
    }

    if (lowerPtr_)
    {
        gpuField<LUType>& lower = *lowerPtr_;

        const labelgpuList& u = lduAddr().upperAddr();

        thrust::transform
        (
            lower.begin(),
            lower.end(),
            thrust::make_permutation_iterator(sf.begin(),u.begin()),
            lower.begin(),
            multiplyOperatorFunctor<LUType,LUType,LUType>()
        );
    }

    FatalErrorIn
    (
        "LduMatrix<Type, DType, LUType>::operator*=(const scalarField& sf)"
    )   << "Scaling a matrix by scalarField is not currently supported\n"
           "because scaling interfacesUpper_ and interfacesLower_ "
           "require special transfers"
        << abort(FatalError);

    //interfacesUpper_ *= ;
    //interfacesLower_ *= sf;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (sourcePtr_)
    {
        *sourcePtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }

    interfacesUpper_ *= s;
    interfacesLower_ *= s;
}


// ************************************************************************* //
