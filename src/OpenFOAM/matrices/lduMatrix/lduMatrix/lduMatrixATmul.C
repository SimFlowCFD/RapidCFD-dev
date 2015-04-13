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

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "textures.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{                
        
#define MAX_NEI_SIZE 3

struct matrixMultiplyFunctor
{
    textures<scalar> psi;
    const scalar* lower;
    const scalar* upper;
    const label* own;
    const label* nei;

    matrixMultiplyFunctor
    (
        textures<scalar> _psi, 
        const scalar* _lower,
        const scalar* _upper,
        const label* _own,
        const label* _nei
    ):
        psi(_psi),
        lower(_lower),
        upper(_upper),
        own(_own),
        nei(_nei)
    {}

    __device__
    scalar operator()(const scalar& d,const thrust::tuple<label,label,label,label>& t)
    {
        scalar out = d;
        scalar tmpSum[2*MAX_NEI_SIZE] = {};
        scalar nExtra = 0;
            
        label oStart = thrust::get<0>(t);
        label oSize = thrust::get<1>(t) - oStart;
            
        label nStart = thrust::get<2>(t);
        label nSize = thrust::get<3>(t) - nStart;

        for(label i = 0; i<MAX_NEI_SIZE; i++)
        {
            if(i<oSize)
            {
                label face = oStart + i;

                tmpSum[i] = upper[face]*psi[nei[face]]; 
            }
        }

        for(label i = 0; i<MAX_NEI_SIZE; i++)
        {
            if(i<nSize)
            {
                 label face = nStart + i;
                   
                 tmpSum[i+MAX_NEI_SIZE] = lower[face]*psi[own[face]];
            }
        }

        for(label i = 0; i<2*MAX_NEI_SIZE; i++)
        {
            out+= tmpSum[i]; 
        }
           
        for(label i = MAX_NEI_SIZE; i<oSize; i++)
        {
            label face = oStart + i;
                
            out += upper[face]*psi[nei[face]];
        }
            
            
        for(label i = MAX_NEI_SIZE; i<nSize; i++)
        {
            label face = nStart + i;

            nExtra += lower[face]*psi[own[face]]; 
        }  
            
        return out + nExtra;
    }
};
    
#undef MAX_NEI_SIZE

inline void callMultiply
(
    scalargpuField& Apsi,
    const scalargpuField& psi,

    const labelgpuList& l,
    const labelgpuList& u,

    const labelgpuList& ownStart,
    const labelgpuList& losortStart,

    const scalargpuField& Lower,
    const scalargpuField& Upper,
    const scalargpuField& Diag
)
{
    textures<scalar> psiTex(psi);

    thrust::transform
    (
        thrust::make_transform_iterator
        (
            thrust::make_zip_iterator(thrust::make_tuple
            (
                Diag.begin(),
                psi.begin()
            )),
            lduMatrixDiagonalFunctor()
        ),
        thrust::make_transform_iterator
        (
            thrust::make_zip_iterator(thrust::make_tuple
            (
                Diag.end(),
                psi.end()
            )),
            lduMatrixDiagonalFunctor()
        ),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            ownStart.begin(),
            ownStart.begin()+1,
            losortStart.begin(),
            losortStart.begin()+1
        )),
        Apsi.begin(),
        matrixMultiplyFunctor
        (
            psiTex,
            Lower.data(),
            Upper.data(),
            l.data(),
            u.data()
        )
    );

    psiTex.destroy();
}


}

void Foam::lduMatrix::Amul
(
    scalargpuField& Apsi,
    const tmp<scalargpuField>& tpsi,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const labelgpuList& l = lduAddr().ownerSortAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const scalargpuField& Lower = lowerSort();
    const scalargpuField& Upper = upper();
    const scalargpuField& Diag = diag();

    const scalargpuField& psi = tpsi();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    callMultiply
    (
        Apsi,
        psi,
        l,
        u,
        ownStart,
        losortStart,
        Lower,
        Upper,
        Diag
    );

    updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    tpsi.clear();
}


void Foam::lduMatrix::Tmul
(
    scalargpuField& Tpsi,
    const tmp<scalargpuField>& tpsi,
    const FieldField<gpuField, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const labelgpuList& l = lduAddr().ownerSortAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();

    const scalargpuField& Lower = lower();
    const scalargpuField& Upper = upperSort();
    const scalargpuField& Diag = diag();

    const scalargpuField& psi = tpsi();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    callMultiply
    (
        Tpsi,
        psi,
        l,
        u,
        ownStart,
        losortStart,
        Upper,
        Lower,
        Diag
    );

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    tpsi.clear();
}


void Foam::lduMatrix::sumA
(
    scalargpuField& sumA,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    matrixOperation
    (
        diag().begin(),
        sumA,
        lduAddr(),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            upper().data(),
            unityOp<scalar>()
        ),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            lowerSort().data(),
            unityOp<scalar>()
        )
    );


    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces, patchI)
    {
        if (interfaces.set(patchI))
        {
            const scalargpuField& pCoeffs = interfaceBouCoeffs[patchI];
                                               
            matrixPatchOperation
            (
                patchI,
                sumA,
                lduAddr(),
                matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >
                (
                    pCoeffs.data(),
                    negateUnaryOperatorFunctor<scalar,scalar>()
                )
            );
        }
    }
}


void Foam::lduMatrix::residual
(
    scalargpuField& rA,
    const scalargpuField& psi,
    const scalargpuField& source,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const labelgpuList& l = lduAddr().ownerSortAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const scalargpuField& Lower = lowerSort();
    const scalargpuField& Upper = upper();
    const scalargpuField& Diag = diag();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<gpuField, scalar> mBouCoeffs(interfaceBouCoeffs.size());

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces.set(patchi))
        {
            mBouCoeffs.set(patchi, -interfaceBouCoeffs[patchi]);
        }
    }

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );
								   
    matrixOperation
    (
        thrust::make_transform_iterator
        (
            thrust::make_zip_iterator(thrust::make_tuple
            ( 
                 source.begin(),
                 Diag.begin(),
                 psi.begin() 
            )), 
            lduMatrixDiagonalResidualFunctor() 
        ),
        rA,
        lduAddr(),
        matrixCoeffsMultiplyFunctor<scalar,scalar,negateUnaryOperatorFunctor<scalar,scalar> >
        (
            psi.data(),
            Upper.data(),
            u.data(),
            negateUnaryOperatorFunctor<scalar,scalar>()
        ),
        matrixCoeffsMultiplyFunctor<scalar,scalar,negateUnaryOperatorFunctor<scalar,scalar> >
        (
            psi.data(),
            Lower.data(),
            l.data(),
            negateUnaryOperatorFunctor<scalar,scalar>()
        )
    );                                       

    // Update interface interfaces
    updateMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );
}


Foam::tmp<Foam::scalargpuField> Foam::lduMatrix::residual
(
    const scalargpuField& psi,
    const scalargpuField& source,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    tmp<scalargpuField> trA(new scalargpuField(psi.size()));
    residual(trA(), psi, source, interfaceBouCoeffs, interfaces, cmpt);
    return trA;
}


Foam::tmp<Foam::scalargpuField > Foam::lduMatrix::H1() const
{
    tmp<scalargpuField > tH1
    (
        new scalargpuField(lduAddr().size(), 0.0)
    );

    if (lowerPtr_ || upperPtr_)
    {
        scalargpuField& H_ = tH1();

        const scalargpuField& Lower = lowerSort();
        const scalargpuField& Upper = upper();
        
        matrixOperation
        (
            H_.begin(),
            H_,
            lduAddr(),
            matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >
            (
                Upper.data(),
                negateUnaryOperatorFunctor<scalar,scalar>()
            ),
            matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >
            (
                Lower.data(),
                negateUnaryOperatorFunctor<scalar,scalar>()
            )
        );

    }

    return tH1;
}


// ************************************************************************* //
