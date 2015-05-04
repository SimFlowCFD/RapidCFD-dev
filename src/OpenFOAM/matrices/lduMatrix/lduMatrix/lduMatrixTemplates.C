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

Description
    lduMatrix member H operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "lduMatrixFunctors.H"
#include "lduAddressingFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
struct lduMatrixfaceHFunctor
{
    template<class Tuple>
    __HOST____DEVICE__
    Type operator()(const Tuple& t)
    {
        return thrust::get<0>(t)*thrust::get<1>(t) 
               - thrust::get<2>(t)*thrust::get<3>(t);
    }
};

}

template<class Type>
void Foam::lduMatrix::H(Foam::gpuField<Type>& Hpsi,const gpuField<Type>& psi) const
{
    Hpsi = pTraits<Type>::zero;

    if (lowerPtr_ || upperPtr_)
    {
        const scalargpuField& Lower = this->lower();
        const scalargpuField& Upper = this->upper();

        const labelgpuList& l = lduAddr().lowerAddr();
        const labelgpuList& u = lduAddr().upperAddr();
        
        matrixOperation
        (
            Hpsi.begin(),
            Hpsi,
            lduAddr(),
            matrixCoeffsMultiplyFunctor<Type,scalar,negateUnaryOperatorFunctor<Type,Type> >
            (
                psi.data(),
                Upper.data(),
                u.data(),
                negateUnaryOperatorFunctor<Type,Type>()
            ),
            matrixCoeffsMultiplyFunctor<Type,scalar,negateUnaryOperatorFunctor<Type,Type> >
            (
                psi.data(),
                Lower.data(),
                l.data(),
                negateUnaryOperatorFunctor<Type,Type>()
            )
        );                                        
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::lduMatrix::H(const gpuField<Type>& psi) const
{
    tmp<gpuField<Type> > tHpsi
    (
        new gpuField<Type>(lduAddr().size(), pTraits<Type>::zero)
    );

    H(tHpsi(),psi);

    return tHpsi;
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::lduMatrix::H(const tmp<gpuField<Type> >& tpsi) const
{
    tmp<gpuField<Type> > tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}

template<class Type>
void Foam::lduMatrix::faceH(Foam::gpuField<Type>& faceHpsi,const gpuField<Type>& psi) const
{
    if (lowerPtr_ || upperPtr_)
    {
        const scalargpuField& Lower = const_cast<const lduMatrix&>(*this).lower();
        const scalargpuField& Upper = const_cast<const lduMatrix&>(*this).upper();

        const labelgpuList& l = lduAddr().lowerAddr();
        const labelgpuList& u = lduAddr().upperAddr();

        thrust::transform
        ( 
            thrust::make_zip_iterator(thrust::make_tuple
            (
                Upper.begin(),  
                thrust::make_permutation_iterator(psi.begin(), u.begin()),
                Lower.begin(),
                thrust::make_permutation_iterator(psi.begin(), l.begin())
            )),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                Upper.end() ,  
                thrust::make_permutation_iterator(psi.begin(), u.end()),
                Lower.end(),
                thrust::make_permutation_iterator(psi.begin(), l.end())
            )),
            faceHpsi.begin(),
            lduMatrixfaceHFunctor<Type>()
        );
    }
    else
    {
        FatalErrorIn("lduMatrix::faceH(const Field<Type>& psi) const")
            << "Cannot calculate faceH"
               " the matrix does not have any off-diagonal coefficients."
            << exit(FatalError);
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::lduMatrix::faceH(const gpuField<Type>& psi) const
{
    tmp<gpuField<Type> > tfaceHpsi(new gpuField<Type> (lower().size()));
    gpuField<Type> & faceHpsi = tfaceHpsi();

    faceH(faceHpsi,psi);

    return tfaceHpsi;
}


template<class Type>
Foam::tmp<Foam::gpuField<Type> >
Foam::lduMatrix::faceH(const tmp<gpuField<Type> >& tpsi) const
{
    tmp<gpuField<Type> > tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// ************************************************************************* //
