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
#include "Textures.H"
#include "lduMatrixSolutionCache.H"

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<bool fast,int nUnroll>
struct matrixMultiplyFunctor
{
    const textures<scalar> psi;
    const scalar * diag;
    const scalar * lower;
    const scalar * upper;
    const label * own;
    const label * nei;
    const label * ownStart;
    const label * losortStart;
    const label * losort;

    matrixMultiplyFunctor
    (
        const textures<scalar> _psi,
        const scalar * _diag,
        const scalar * _lower,
        const scalar * _upper,
        const label * _own,
        const label * _nei,
        const label * _ownStart,
        const label * _losortStart,
        const label * _losort
    ):
        psi(_psi),
        diag(_diag),
        lower(_lower),
        upper(_upper),
        own(_own),
        nei(_nei),
        ownStart(_ownStart),
        losortStart(_losortStart),
        losort(_losort)
    {}

    __device__
    scalar operator()(const label& id) const
    {
        scalar tmpSum[2*nUnroll] = {};
        scalar nExtra = 0;

        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;

        label nStart = losortStart[id];
        label nSize = losortStart[id+1] - nStart;

        scalar out = diag[id]*psi[id];

        for(label i = 0; i<nUnroll; i++)
        {
            if(i<oSize)
            {
                label face = oStart + i;

                tmpSum[i] = upper[face]*psi[nei[face]];
            }
        }

        for(label i = 0; i<nUnroll; i++)
        {
            if(i<nSize)
            {
                 label face = nStart + i;
                 if( ! fast)
                     face = losort[face];

                 tmpSum[i+nUnroll] = lower[face]*psi[own[face]];
            }
        }

        #pragma unroll
        for(label i = 0; i<2*nUnroll; i++)
        {
            out+= tmpSum[i];
        }

        for(label i = nUnroll; i<oSize; i++)
        {
            label face = oStart + i;

            out += upper[face]*psi[nei[face]];
        }

        for(label i = nUnroll; i<nSize; i++)
        {
            label face = nStart + i;
            if( ! fast)
                face = losort[face];

            nExtra += lower[face]*psi[own[face]];
        }

        return out + nExtra;
    }
};

template<bool fast>
inline void callMultiply
(
    scalargpuField& Apsi,
    const scalargpuField& psi,

    const labelgpuList& l,
    const labelgpuList& u,

    const labelgpuList& ownStart,
    const labelgpuList& losortStart,
    const labelgpuList& losort,

    const scalargpuField& Lower,
    const scalargpuField& Upper,
    const scalargpuField& Diag
)
{
    textureBind<scalar> psiTex(psi);

    thrust::transform
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+psi.size(),
        Apsi.begin(),
        matrixMultiplyFunctor<fast,3>
        (
            psiTex(),
            Diag.data(),
            Lower.data(),
            Upper.data(),
            l.data(),
            u.data(),
            ownStart.data(),
            losortStart.data(),
            losort.data()
        )
    );
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
    bool fastPath = lduMatrixSolutionCache::favourSpeed >= 2 ||
                    (lduMatrixSolutionCache::favourSpeed && ( coarsestLevel() || ! level()));

    const labelgpuList& l = fastPath? lduAddr().ownerSortAddr(): lduAddr().lowerAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();
    const labelgpuList& losort = lduAddr().losortAddr();

    const scalargpuField& Lower = fastPath? lowerSort(): lower();
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

    if(fastPath)
    {
        callMultiply<true>
        (
            Apsi,
            psi,
            l,
            u,
            ownStart,
            losortStart,
            losort,
            Lower,
            Upper,
            Diag
        );
    }
    else
    {
        callMultiply<false>
        (
            Apsi,
            psi,
            l,
            u,
            ownStart,
            losortStart,
            losort,
            Lower,
            Upper,
            Diag
        );
    }

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
    bool fastPath = lduMatrixSolutionCache::favourSpeed;

    const labelgpuList& l = fastPath? lduAddr().ownerSortAddr(): lduAddr().lowerAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const labelgpuList& ownStart = lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = lduAddr().losortStartAddr();
    const labelgpuList& losort = lduAddr().losortAddr();

    const scalargpuField& Lower = lower();
    const scalargpuField& Upper = fastPath? upperSort(): upper();
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

    if(fastPath)
    {
        callMultiply<true>
        (
            Tpsi,
            psi,
            l,
            u,
            ownStart,
            losortStart,
            losort,
            Upper,
            Lower,
            Diag
        );
    }
    else
    {
        callMultiply<false>
        (
            Tpsi,
            psi,
            l,
            u,
            ownStart,
            losortStart,
            losort,
            Upper,
            Lower,
            Diag
        );
    }

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
    const scalargpuField& Lower = lower();
    const scalargpuField& Upper = upper();
    const scalargpuField& Diag = diag();

    matrixOperation
    (
        Diag.begin(),
        sumA,
        lduAddr(),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            Upper.data(),
            unityOp<scalar>()
        ),
        matrixCoeffsFunctor<scalar,unityOp<scalar> >
        (
            Lower.data(),
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

#define CALL_RESIDUAL_FUNCTION(functionName)                                                \
functionName                                                                                \
(                                                                                           \
    thrust::make_transform_iterator                                                         \
    (                                                                                       \
        thrust::make_zip_iterator(thrust::make_tuple                                        \
        (                                                                                   \
            source.begin(),                                                                 \
            Diag.begin(),                                                                   \
            psi.begin()                                                                     \
        )),                                                                                 \
        lduMatrixDiagonalResidualFunctor()                                                  \
    ),                                                                                      \
    rA,                                                                                     \
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
    bool fastPath = lduMatrixSolutionCache::favourSpeed;

    const labelgpuList& l = fastPath? lduAddr().ownerSortAddr(): lduAddr().lowerAddr();
    const labelgpuList& u = lduAddr().upperAddr();

    const scalargpuField& Lower = fastPath? lowerSort(): lower();
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

    if(fastPath)
    {
        CALL_RESIDUAL_FUNCTION(matrixFastOperation);
    }
    else
    {
        CALL_RESIDUAL_FUNCTION(matrixOperation);
    }

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

#undef CALL_RESIDUAL_FUNCTION


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

#define CALL_H_FUNCTION(functionName)                                                \
functionName                                                                         \
(                                                                                    \
    H_.begin(),                                                                      \
    H_,                                                                              \
    lduAddr(),                                                                       \
    matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >           \
    (                                                                                \
                Upper.data(),                                                        \
                negateUnaryOperatorFunctor<scalar,scalar>()                          \
    ),                                                                               \
    matrixCoeffsFunctor<scalar,negateUnaryOperatorFunctor<scalar,scalar> >           \
    (                                                                                \
        Lower.data(),                                                                \
        negateUnaryOperatorFunctor<scalar,scalar>()                                  \
    )                                                                                \
)

void Foam::lduMatrix::H1(scalargpuField& H_) const
{
    H_ = 0.0;

    if (lowerPtr_ || upperPtr_)
    {
        bool fastPath = lduMatrixSolutionCache::favourSpeed;

        const scalargpuField& Lower = fastPath?lowerSort():lower();
        const scalargpuField& Upper = upper();

        if(fastPath)
        {
            CALL_H_FUNCTION(matrixFastOperation);
        }
        else
        {
            CALL_H_FUNCTION(matrixOperation);
        }

    }
}

Foam::tmp<Foam::scalargpuField > Foam::lduMatrix::H1() const
{
    tmp<scalargpuField > tH1
    (
        new scalargpuField(lduAddr().size(), 0.0)
    );

    H1(tH1());

    return tH1;
}

#undef CALL_H_FUNCTION


// ************************************************************************* //
