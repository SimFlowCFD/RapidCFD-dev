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

#include "LduMatrix.H"
#include "LduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class LUType>
class Amultiplier
:
    public LduInterfaceField<Type>::Amultiplier
{
    const gpuField<LUType>& A_;

public:

    Amultiplier(const gpuField<LUType>& A)
    :
        A_(A)
    {}

    virtual ~Amultiplier()
    {}

    virtual void addAmul(gpuField<Type>& Apsi, const gpuField<Type>& psi) const
    {
        Apsi += A_*psi;
    }
};

template<class Type,class DType>
struct LduMatrixDiagonalResidualFunctor : public std::unary_function<thrust::tuple<Type,DType,Type>,Type>
{
    __HOST____DEVICE__
    Type operator()(const thrust::tuple<Type,DType,Type>& c)
    {
        return thrust::get<0>(c) -       //source 
               dot(thrust::get<1>(c),    //diagonal
                   thrust::get<2>(c));   //psi
    }
};

    
template<class Type,class LUType>
struct LduMatrixPatchSubtractFunctor : public std::binary_function<label,Type,Type>
{
    const Type one;
    const LUType* coeff;
    const label* neiStart;
    const label* losort;

    LduMatrixPatchSubtractFunctor
    (
        const Type _one,
        const LUType* _coeff,
        const label* _neiStart,
        const label* _losort
    ):
        one(_one),
        coeff(_coeff),
        neiStart(_neiStart),
        losort(_losort)
    {}

    __HOST____DEVICE__
    Type operator()(const label& id,const Type& s)
    {
        Type out = s;

        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            out -= dot(coeff[face], one);
        }

        return out;
    }
};
    
template<class Type, class DType, class LUType,bool addOffDiagonal, bool normalMult>
struct LduMatrixMultiplyFunctor : public std::binary_function<Type,label,Type>
{
    const Type* psi;
    const LUType* lower;
    const LUType* upper;
    const label* ownStart;
    const label* neiStart;
    const label* own;
    const label* nei;
    const label* losort;

    LduMatrixMultiplyFunctor
    (
        const Type* _psi, 
        const LUType* _lower,
        const LUType* _upper,
        const label* _ownStart,
        const label* _neiStart,
        const label* _own,
        const label* _nei,
        const label* _losort
    ):
        psi(_psi),
        lower(_lower),
        upper(_upper),
        ownStart(_ownStart),
        neiStart(_neiStart),
        own(_own),
        nei(_nei),
        losort(_losort)
    {}

    __HOST____DEVICE__
    Type operator()(const Type& d,const label& id)
    {
        Type out = d;
        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;
            
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<oSize; i++)
        {
            label face = oStart + i;
            if(addOffDiagonal)
            {
                if(normalMult)
                    out += dot(upper[face],psi[nei[face]]); 
                else
                    out += dot(lower[face],psi[nei[face]]); 
            }
            else
            {
                if(normalMult)
                    out -= dot(upper[face],psi[nei[face]]); 
                else
                    out -= dot(lower[face],psi[nei[face]]); 
            }
        }


        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            if(addOffDiagonal)
            {
                if(normalMult)
                    out += dot(lower[face],psi[own[face]]); 
                else
                    out += dot(upper[face],psi[own[face]]);
            }
            else
            {
                if(normalMult)
                    out -= dot(lower[face],psi[own[face]]); 
                else
                    out -= dot(upper[face],psi[own[face]]); 
            }
        }

        return out;
    }
};

template<class Type, class DType, class LUType,bool addOffDiagonal>
struct LduMatrixSumFunctor : public std::binary_function<Type,label,Type>
{
    const Type one;
    const LUType* lower;
    const LUType* upper;
    const label* ownStart;
    const label* neiStart;
    const label* own;
    const label* nei;
    const label* losort;

    LduMatrixSumFunctor
    (
        const Type _one,
        const LUType* _lower,
        const LUType* _upper,
        const label* _ownStart,
        const label* _neiStart,
        const label* _own,
        const label* _nei,
        const label* _losort
    ):
        one(_one),
        lower(_lower),
        upper(_upper),
        ownStart(_ownStart),
        neiStart(_neiStart),
        own(_own),
        nei(_nei),
        losort(_losort)
    {}

    __HOST____DEVICE__
    Type operator()(const Type& diag,const label& id)
    {
        Type out = diag;
        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;

        for(label i = 0; i<oSize; i++)
        {
            label face = oStart + i;
            if(addOffDiagonal)
            {
                out += dot(upper[face],one);
            }
            else
            {
                out -= dot(upper[face],one);
            }
        }

        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];
            if(addOffDiagonal)
            {
                out += dot(lower[face],one);
            }
            else
            {
                out -= dot(lower[face],one); 
            }
        }

        return out;
    }
};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::Amul
(
    gpuField<Type>& Apsi,
    const tmp<gpuField<Type> >& tpsi
) const
{
    
    const gpuField<LUType>& Upper = upper();
    const gpuField<LUType>& Lower = lower();
    const gpuField<DType>& Diag = diag();
    const gpuField<Type>& psi = tpsi();
    
    const labelgpuList& losort = lduAddr().losortAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const labelgpuList& u = lduAddr().upperAddr();
    const labelgpuList& l = lduAddr().lowerAddr();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfacesUpper_,
        psi,
        Apsi
    );

    thrust::transform
    (
        Diag.begin(),
        Diag.end(),
        psi.begin(),
        Apsi.begin(),
        dotBinaryFunctionFunctor<DType,Type,Type>()
    );

   thrust::transform
   (
        Apsi.begin(),
        Apsi.end(),
        thrust::make_counting_iterator(0),
        Apsi.begin(),
        LduMatrixMultiplyFunctor<Type,DType,LUType,true,true>
        (
            psi.data(),
            Lower.data(),
            Upper.data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfacesUpper_,
        psi,
        Apsi
    );

    tpsi.clear();
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::Tmul
(
    gpuField<Type>& Tpsi,
    const tmp<gpuField<Type> >& tpsi
) const
{
    const gpuField<LUType>& Upper = upper();
    const gpuField<LUType>& Lower = lower();
    const gpuField<DType>& Diag = diag();
    const gpuField<Type>& psi = tpsi();
    
    const labelgpuList& losort = lduAddr().losortAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const labelgpuList& u = lduAddr().upperAddr();
    const labelgpuList& l = lduAddr().lowerAddr();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfacesLower_,
        psi,
        Tpsi
    );

    thrust::transform
    (
        Diag.begin(),
        Diag.end(),
        psi.begin(),
        Tpsi.begin(),
        dotBinaryFunctionFunctor<DType,Type,Type>()
    );
                       
   thrust::transform
   (
        Tpsi.begin(),
        Tpsi.end(),
        thrust::make_counting_iterator(0),
        Tpsi.begin(),
        LduMatrixMultiplyFunctor<Type,DType,LUType,true,false>
        (
            psi.data(),
            Lower.data(),
            Upper.data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfacesLower_,
        psi,
        Tpsi
    );

    tpsi.clear();
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumA
(
    gpuField<Type>& sumA
) const
{
    const gpuField<LUType>& Upper = upper();
    const gpuField<LUType>& Lower = lower();
    const gpuField<DType>& Diag = diag();
    
    const labelgpuList& losort = lduAddr().losortAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const labelgpuList& u = lduAddr().upperAddr();
    const labelgpuList& l = lduAddr().lowerAddr();

    thrust::transform
    (
        Diag.begin(),
        Diag.end(),
        sumA.begin(),
        dotBinaryFunctionFSFunctor<DType,Type,Type>(pTraits<Type>::one)
    );
             
    thrust::transform
    (
        sumA.begin(),
        sumA.end(),
        thrust::make_counting_iterator(0),
        sumA.begin(),
        LduMatrixSumFunctor<Type,DType,LUType,true>
        (
            pTraits<Type>::one,
            Lower.data(),
            Upper.data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces_, patchI)
    {
        if (interfaces_.set(patchI))
        {
            const labelgpuList& pcells = lduAddr().patchSortCells(patchI);
            const gpuField<LUType>& pCoeffs = interfacesUpper_[patchI];
            
            const labelgpuList& losort = lduAddr().patchSortAddr(patchI);
            const labelgpuList& losortStart = lduAddr().patchSortStartAddr(patchI);
  
            thrust::transform
            (
                thrust::make_counting_iterator(0),
                thrust::make_counting_iterator(0)+pcells.size(),
                thrust::make_permutation_iterator(sumA.begin(),pcells.begin()),
                thrust::make_permutation_iterator(sumA.begin(),pcells.begin()),
                LduMatrixPatchSubtractFunctor<Type,LUType>
                (
                    pTraits<Type>::one,
                    pCoeffs.data(),
                    losortStart.data(),
                    losort.data()
                )
            );
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::residual
(
    gpuField<Type>& rA,
    const gpuField<Type>& psi
) const
{
    const gpuField<LUType>& Upper = upper();
    const gpuField<LUType>& Lower = lower();
    const gpuField<DType>& Diag = diag();
    const gpuField<Type>& Source = source();
    
    const labelgpuList& losort = lduAddr().losortAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const labelgpuList& u = lduAddr().upperAddr();
    const labelgpuList& l = lduAddr().lowerAddr();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update to add the contibution to the r.h.s.

    FieldField<gpuField, LUType> mBouCoeffs(interfacesUpper_.size());

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs.set(patchi, -interfacesUpper_[patchi]);
        }
    }

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        mBouCoeffs,
        psi,
        rA
    );
        
    thrust::transform
    (
        thrust::make_transform_iterator
        (
            thrust::make_zip_iterator( thrust::make_tuple
            ( 
                Source.begin(),
                Diag.begin(),
                psi.begin()
            )), 
            LduMatrixDiagonalResidualFunctor<Type,DType>() 
        ),
        thrust::make_transform_iterator
        (
            thrust::make_zip_iterator( thrust::make_tuple
            ( 
                Source.end(),
                Diag.end(),
                psi.end() 
            )), 
           LduMatrixDiagonalResidualFunctor<Type,DType>() 
        ),
        thrust::make_counting_iterator(0),
        rA.begin(),
        LduMatrixMultiplyFunctor<Type,DType,LUType,false,true>
        (
            psi.data(),
            Lower.data(),
            Upper.data(),
            ownStart.data(),
            losortStart.data(),
            l.data(),
            u.data(),
            losort.data()
        )
    );

    // Update interface interfaces
    updateMatrixInterfaces
    (
        mBouCoeffs,
        psi,
        rA
    );
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::gpuField<Type> > Foam::LduMatrix<Type, DType, LUType>::residual
(
    const gpuField<Type>& psi
) const
{
    tmp<gpuField<Type> > trA(new gpuField<Type>(psi.size()));
    residual(trA(), psi);
    return trA;
}


// ************************************************************************* //
