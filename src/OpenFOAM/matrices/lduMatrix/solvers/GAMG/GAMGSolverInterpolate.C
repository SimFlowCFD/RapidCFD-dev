/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "GAMGSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

struct GAMGinterpolatePsiFunctor
{
    __HOST____DEVICE__
    scalar operator()(const scalar& Apsi, const scalar& diag)
    {
        return - Apsi/diag;
    }
};

}


void Foam::GAMGSolver::interpolate
(
    scalargpuField& psi,
    scalargpuField& Apsi,
    const lduMatrix& m,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    const labelgpuList& l = m.lduAddr().ownerSortAddr();
    const labelgpuList& u = m.lduAddr().upperAddr();

    const scalargpuField& Lower = m.lowerSort();
    const scalargpuField& Upper = m.upper();
    const scalargpuField& Diag = m.diag();

    m.initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    matrixFastOperation
    (
        thrust::make_constant_iterator(scalar(0.0)),
        Apsi,
        m.lduAddr(),
        matrixCoeffsMultiplyFunctor<scalar,scalar,thrust::identity<scalar> >
        (
            psi.data(),
            Upper.data(),
            u.data(),
            thrust::identity<scalar>()
        ),
        matrixCoeffsMultiplyFunctor<scalar,scalar,thrust::identity<scalar> >
        (
            psi.data(),
            Lower.data(),
            l.data(),
            thrust::identity<scalar>()
        )
    );  

    m.updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    thrust::transform
    (
        Apsi.begin(),
        Apsi.end(),
        Diag.begin(),
        psi.begin(),
        GAMGinterpolatePsiFunctor()
    );
}

namespace Foam
{

struct GAMGinterpolateCorrCdiagCFunctor
{
    const scalar* diag;
    const scalar* psi;
    const label* sort;

    GAMGinterpolateCorrCdiagCFunctor
    (
        const scalar* _diag,
        const scalar* _psi,
        const label* _sort
    ):
        diag(_diag),
        psi(_psi),
        sort(_sort)
    {}

    __HOST____DEVICE__
    thrust::tuple<scalar,scalar> operator()
    (
        const label& start, 
        const label& end
    )
    {
        scalar corrC = 0;
        scalar diagC = 0;

        for(label i = start; i<end; i++)
        {
            label celli = sort[i];

            corrC += diag[celli]*psi[celli];
            diagC += diag[celli];
        }

        return thrust::make_tuple(corrC,diagC);
    }

};

struct GAMGinterpolateCorrCFunctor
{
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& psiC, const Tuple& t)
    {
        return psiC - thrust::get<0>(t)/thrust::get<1>(t);
    }
}; 
   
}

void Foam::GAMGSolver::interpolate
(
    scalargpuField& psi,
    scalargpuField& Apsi,
    const lduMatrix& m,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const labelgpuList& restrictSortAddressing,
    const labelgpuList& restrictTargetAddressing,
    const labelgpuList& restrictTargetStartAddressing,
    const scalargpuField& psiC,
    const direction cmpt
) const
{

    notImplemented("GAMGSolver::interpolate()");

    interpolate
    (
        psi,
        Apsi,
        m,
        interfaceBouCoeffs,
        interfaces,
        cmpt
    );

    const label nCCells = psiC.size();
    scalargpuField corrC(nCCells, 0);
    scalargpuField diagC(nCCells, 0);

    thrust::transform
    (
        restrictTargetStartAddressing.begin(),
        restrictTargetStartAddressing.end()-1,
        restrictTargetStartAddressing.begin()+1,
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator
            (
                corrC.begin(),
                restrictTargetAddressing.begin()
            ),
            thrust::make_permutation_iterator
            (
                diagC.begin(),
                restrictTargetAddressing.begin()
            )
        )),
        GAMGinterpolateCorrCdiagCFunctor
        (
            m.diag().data(),
            psi.data(),
            restrictSortAddressing.data()
        )
    );
   
    thrust::transform
    (
        psiC.begin(),
        psiC.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            corrC.begin(),
            diagC.begin()
        )),
        corrC.begin(),
        GAMGinterpolateCorrCFunctor()
    );
    //TODO
    /*
    thrust::transform
    (
        psi.begin(),
        psi.end(),
        thrust::make_permutation_iterator
        (
            corrC.begin(),
            restrictAddressing.begin()
        ),
        psi.begin(),
        thrust::plus<scalar>()
    );
*/
}


// ************************************************************************* //
