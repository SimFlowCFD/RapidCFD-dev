/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "CMULES.H"
#include "fvcSurfaceIntegrate.H"
#include "slicedSurfaceFields.H"
#include "wedgeFvPatch.H"
#include "syncTools.H"
#include "MULESFunctors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::correct
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su
)
{
    Info<< "MULES: Correcting " << psi.name() << endl;

    const fvMesh& mesh = psi.mesh();

    scalargpuField psiIf(psi.size(), 0);
    fvc::surfaceIntegrate(psiIf, phiCorr);

    if (mesh.moving())
    {
        psi.internalField() =
        (
            rho.getField()*psi.internalField()*rDeltaT
          + Su.getField()
          - psiIf
        )/(rho.getField()*rDeltaT - Sp.getField());
    }
    else
    {
        psi.internalField() =
        (
            rho.getField()*psi.internalField()*rDeltaT
          + Su.getField()
          - psiIf
        )/(rho.getField()*rDeltaT - Sp.getField());
    }

    psi.correctBoundaryConditions();
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const fvMesh& mesh = psi.mesh();
    const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

    const dictionary& MULEScontrols = mesh.solverDict(psi.name());

    label nLimiterIter
    (
        readLabel(MULEScontrols.lookup("nLimiterIter"))
    );

    limitCorr
    (
        rDeltaT,
        rho,
        psi,
        phi,
        phiCorr,
        Sp, Su,
        psiMax, psiMin,
        nLimiterIter
    );

    correct(rDeltaT, rho, psi, phi, phiCorr, Sp, Su);
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::LTScorrect
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    const volScalarField& rDeltaT =
        mesh.objectRegistry::lookupObject<volScalarField>("rSubDeltaT");

    const dictionary& MULEScontrols = mesh.solverDict(psi.name());

    label nLimiterIter
    (
        readLabel(MULEScontrols.lookup("nLimiterIter"))
    );

    limitCorr
    (
        rDeltaT,
        rho,
        psi,
        phi,
        phiCorr,
        Sp, Su,
        psiMax, psiMin,
        nLimiterIter
    );
    correct(rDeltaT, rho, psi, phi, phiCorr, Sp, Su);
}



namespace Foam
{

struct limiterCMULESFunctor
{
    const label* own;
    const label* nei;
    const label* ownStart;
    const label* neiStart;
    const label* losort;

    const scalar* psiIf;
    const scalar* phiCorrIf;

    scalar* psiMaxn;
    scalar* psiMinn;
    scalar* sumPhip;
    scalar* mSumPhim;

    limiterCMULESFunctor
    (
        const label* _own,
        const label* _nei,
        const label* _ownStart,
        const label* _neiStart,
        const label* _losort,

        const scalar* _psiIf,
        const scalar* _phiCorrIf,

        scalar* _psiMaxn,
        scalar* _psiMinn,
        scalar* _sumPhip,
        scalar* _mSumPhim
    ):
        own(_own),
        nei(_nei),
        ownStart(_ownStart),
        neiStart(_neiStart),
        losort(_losort),

        psiIf(_psiIf),
        phiCorrIf(_phiCorrIf),

        psiMaxn(_psiMaxn),
        psiMinn(_psiMinn),
        sumPhip(_sumPhip),
        mSumPhim(_mSumPhim)
    {}

    __HOST____DEVICE__
    void operator()(const label& id)
    {
        label oStart = ownStart[id];
        label oSize = ownStart[id+1] - oStart;
		
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        scalar psiMin = psiMinn[id];
        scalar psiMax = psiMaxn[id];
        scalar sumPhipTmp = VSMALL;
        scalar mSumPhimTmp = VSMALL;

        for(label i = 0; i<oSize; i++)
        {
            label face = oStart + i;

            psiMax = max(psiMax,psiIf[nei[face]]);
            psiMin = min(psiMin,psiIf[nei[face]]);

            scalar phiCorrf = phiCorrIf[face];
            if(phiCorrf > 0.0)
            {
                sumPhipTmp += phiCorrf;
            }
            else
            {
                mSumPhimTmp -= phiCorrf;
            }
        }

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];

            psiMax = max(psiMax,psiIf[own[face]]);
            psiMin = min(psiMin,psiIf[own[face]]);

            scalar phiCorrf = phiCorrIf[face];
            if(phiCorrf > 0.0)
            {
                mSumPhimTmp += phiCorrf;
            }
            else
            {
                sumPhipTmp -= phiCorrf;
            }
        }

        psiMaxn[id] = psiMax;
        psiMinn[id] = psiMin;
        sumPhip[id]  = sumPhipTmp;
        mSumPhim[id] = mSumPhimTmp;
    }
};


struct patchMinMaxCMULESFunctor
{
    const label* neiStart;
    const label* losort;
    const label* pCell;

    const scalar* psiPf;
    const scalar* phiCorrPf;

    scalar* psiMaxn;
    scalar* psiMinn;
    scalar* sumPhip;
    scalar* mSumPhim;

    patchMinMaxCMULESFunctor
    (
        const label* _neiStart,
        const label* _losort,
        const label* _pCell,

        const scalar* _psiPf,
        const scalar* _phiCorrPf,

        scalar* _psiMaxn,
        scalar* _psiMinn,
        scalar* _sumPhip,
        scalar* _mSumPhim
    ):
        neiStart(_neiStart),
        losort(_losort),
        pCell(_pCell),

        psiPf(_psiPf),
        phiCorrPf(_phiCorrPf),

        psiMaxn(_psiMaxn),
        psiMinn(_psiMinn),
        sumPhip(_sumPhip),
        mSumPhim(_mSumPhim)
    {}

    __HOST____DEVICE__
    void operator()(const label& id)
    {
        label nStart = neiStart[id];
        label nSize = neiStart[id+1] - nStart;

        label globalId = pCell[id];

        scalar psiMax = psiMaxn[globalId];
        scalar psiMin = psiMinn[globalId];
        scalar sumPhipTmp = sumPhip[globalId];
        scalar mSumPhimTmp = mSumPhim[globalId];

        for(label i = 0; i<nSize; i++)
        {
            label face = losort[nStart + i];

            psiMax = max(psiMax,psiPf[face]);
            psiMin = min(psiMin,psiPf[face]);

            scalar phiCorrf = phiCorrPf[face];

            if (phiCorrf > 0.0)
            {
                sumPhipTmp += phiCorrf;
            }
            else
            {
                mSumPhimTmp -= phiCorrf;
            }
        }

        psiMaxn[globalId] = psiMax;
        psiMinn[globalId] = psiMin;
        sumPhip[globalId] = sumPhipTmp;
        mSumPhim[globalId] = mSumPhimTmp;
    }
};

struct patchLambdaPfCMULESFunctor
{
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& phiCorrfPf, const Tuple& t)
    {
        const scalar& lambdaPf = thrust::get<0>(t);
        const scalar& phiPf = thrust::get<1>(t);

        scalar out = lambdaPf;

        // Limit outlet faces only
        if (phiPf > SMALL*SMALL)
        {
            if (phiCorrfPf > 0.0)
            {
                out = min(lambdaPf, thrust::get<2>(t));
            }
            else
            {
                out = min(lambdaPf, thrust::get<3>(t));
            }
        }
        
        return out;
    }
};

}

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::limiterCorr
(
    scalargpuField& allLambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin,
    const label nLimiterIter
)
{
    const scalargpuField& psiIf = psi;
    const volScalarField::GeometricBoundaryField& psiBf = psi.boundaryField();

    const fvMesh& mesh = psi.mesh();

    const dictionary& MULEScontrols = mesh.solverDict(psi.name());

    scalar extremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>("extremaCoeff", 0.0)
    );

    const labelgpuList& owner = mesh.owner();
    const labelgpuList& neighb = mesh.neighbour();
    const labelgpuList& losort = mesh.lduAddr().losortAddr();

    const labelgpuList& ownStart = mesh.lduAddr().ownerStartAddr();
    const labelgpuList& losortStart = mesh.lduAddr().losortStartAddr();

    tmp<volScalarField::DimensionedInternalField> tVsc = mesh.Vsc();
    const scalargpuField& V = tVsc();

    const surfaceScalarField::GeometricBoundaryField& phiBf =
        phi.boundaryField();

    const scalargpuField& phiCorrIf = phiCorr;
    const surfaceScalarField::GeometricBoundaryField& phiCorrBf =
        phiCorr.boundaryField();

    slicedSurfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless,
        allLambda,
        false   // Use slices for the couples
    );

    scalargpuField& lambdaIf = lambda;
    surfaceScalarField::GeometricBoundaryField& lambdaBf =
        lambda.boundaryField();

    scalargpuField psiMaxn(psiIf.size(), psiMin);
    scalargpuField psiMinn(psiIf.size(), psiMax);

    scalargpuField sumPhip(psiIf.size(), VSMALL);
    scalargpuField mSumPhim(psiIf.size(), VSMALL);

    thrust::for_each
    (
        thrust::make_counting_iterator(0),
        thrust::make_counting_iterator(0)+psiIf.size(),
        limiterCMULESFunctor
        (
            owner.data(),
            neighb.data(),
            ownStart.data(),
            losortStart.data(),
            losort.data(),
            psiIf.data(),
            phiCorrIf.data(),
            psiMaxn.data(),
            psiMinn.data(),
            sumPhip.data(),
            mSumPhim.data()
        )
    );

    forAll(phiCorrBf, patchi)
    {
        const fvPatchScalarField& psiPf = psiBf[patchi];
        const scalargpuField& phiCorrPf = phiCorrBf[patchi];

        const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchi);

        const labelgpuList& losort = mesh.lduAddr().patchSortAddr(patchi);
        const labelgpuList& losortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        thrust::for_each
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+pcells.size(),
            patchMinMaxCMULESFunctor
            (
                losortStart.data(),
                losort.data(),
                pcells.data(),
                psiPf.coupled()?psiPf.patchNeighbourField()().data():psiPf.data(),
                phiCorrPf.data(),
                psiMaxn.data(),
                psiMinn.data(),
                sumPhip.data(),
                mSumPhim.data()
            )
        );
    }

    psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
    psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);

    // scalar smooth = 0.5;
    // psiMaxn = min((1.0 - smooth)*psiIf + smooth*psiMaxn, psiMax);
    // psiMinn = max((1.0 - smooth)*psiIf + smooth*psiMinn, psiMin);

    psiMaxn =
        V
       *(
           (rho.getField()*rDeltaT - Sp.getField())*psiMaxn
         - Su.getField()
         - rho.getField()*psi.internalField()*rDeltaT
        );

    psiMinn =
        V
       *(
           Su.getField()
         - (rho.getField()*rDeltaT - Sp.getField())*psiMinn
         + rho.getField()*psi.internalField()*rDeltaT
        );

    scalargpuField sumlPhip(psiIf.size());
    scalargpuField mSumlPhim(psiIf.size());

    for (int j=0; j<nLimiterIter; j++)
    {
        sumlPhip = 0.0;
        mSumlPhim = 0.0;

        thrust::for_each
        (
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0)+sumlPhip.size(),
            sumlPhiMULESFunctor
            (
                owner.data(),
                neighb.data(),
                ownStart.data(),
                losortStart.data(),
                losort.data(),
                lambdaIf.data(),
                phiCorrIf.data(),
                sumlPhip.data(),
                mSumlPhim.data()
            )
        );

        forAll(lambdaBf, patchi)
        {
            scalargpuField& lambdaPf = lambdaBf[patchi];
            const scalargpuField& phiCorrfPf = phiCorrBf[patchi];

            const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchi);

            const labelgpuList& losort = mesh.lduAddr().patchSortAddr(patchi);
            const labelgpuList& losortStart = mesh.lduAddr().patchSortStartAddr(patchi);

            thrust::for_each
            (
                thrust::make_counting_iterator(0),
                thrust::make_counting_iterator(0)+pcells.size(),
                patchSumlPhiMULESFunctor
                (
                    losortStart.data(),
                    losort.data(),
                    pcells.data(),
                    lambdaPf.data(),
                    phiCorrfPf.data(),
                    sumlPhip.data(),
                    mSumlPhim.data()
                )
            );
        }

        thrust::transform
        (
            sumlPhip.begin(),
            sumlPhip.end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                psiMaxn.begin(),
                mSumPhim.begin()
            )),
            sumlPhip.begin(),
            sumlPhipFinalMULESFunctor<false>()
        );

        thrust::transform
        (
            mSumlPhim.begin(),
            mSumlPhim.end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                psiMinn.begin(),
                sumPhip.begin()
            )),
            mSumlPhim.begin(),
            sumlPhipFinalMULESFunctor<true>()
        );

        const scalargpuField& lambdam = sumlPhip;
        const scalargpuField& lambdap = mSumlPhim;

        thrust::transform
        (
            phiCorrIf.begin(),
            phiCorrIf.end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                lambdaIf.begin(),
                thrust::make_permutation_iterator
                (
                    lambdap.begin(),
                    owner.begin()
                ),
                thrust::make_permutation_iterator
                (
                    lambdam.begin(),
                    neighb.begin()
                ),
                thrust::make_permutation_iterator
                (
                    lambdam.begin(),
                    owner.begin()
                ),
                thrust::make_permutation_iterator
                (
                    lambdap.begin(),
                    neighb.begin()
                )
            )),
            lambdaIf.begin(),
            lambdaIfMULESFunctor()
        );

        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
            const scalargpuField& phiCorrfPf = phiCorrBf[patchi];
            const fvPatchScalarField& psiPf = psiBf[patchi];

            if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
            {
                lambdaPf = 0;
            }
            else if (psiPf.coupled())
            {
                const labelgpuList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                thrust::transform
                (
                    phiCorrfPf.begin(),
                    phiCorrfPf.end(),
                    thrust::make_zip_iterator(thrust::make_tuple
                    (
                        lambdaPf.begin(),
                        thrust::make_permutation_iterator
                        (
                            lambdap.begin(),
                            pFaceCells.begin()
                        ),
                        thrust::make_permutation_iterator
                        (
                            lambdam.begin(),
                            pFaceCells.begin()
                        )
                    )),
                    lambdaPf.begin(),
                    coupledPatchLambdaPfMULESFunctor()
                );
            }
            else
            {
                const labelgpuList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();
                const scalargpuField& phiPf = phiBf[patchi];

                thrust::transform
                (
                    phiCorrfPf.begin(),
                    phiCorrfPf.end(),
                    thrust::make_zip_iterator(thrust::make_tuple
                    (
                        lambdaPf.begin(),
                        phiPf.begin(),
                        thrust::make_permutation_iterator
                        (
                            lambdap.begin(),
                            pFaceCells.begin()
                        ),
                        thrust::make_permutation_iterator
                        (
                            lambdam.begin(),
                            pFaceCells.begin()
                        )
                    )),
                    lambdaPf.begin(),
                    patchLambdaPfCMULESFunctor()
                );
            }
        }

        syncTools::syncFaceList(mesh, allLambda, minOp<scalar>());
    }
}


template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::limitCorr
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin,
    const label nLimiterIter
)
{
    const fvMesh& mesh = psi.mesh();

    scalargpuField allLambda(mesh.nFaces(), 1.0);

    slicedSurfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless,
        allLambda,
        false   // Use slices for the couples
    );

    limiterCorr
    (
        allLambda,
        rDeltaT,
        rho,
        psi,
        phi,
        phiCorr,
        Sp,
        Su,
        psiMax,
        psiMin,
        nLimiterIter
    );

    phiCorr *= lambda;
}


// ************************************************************************* //
