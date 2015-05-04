/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fvScalarMatrix.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMatrixCache.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvMatrix<Foam::scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            scalar delta = diag().get(psi_.mesh().boundary()[patchi].faceCells().get(facei));
			
            internalCoeffs_[patchi].set(facei, internalCoeffs_[patchi].get(facei)+delta);

            boundaryCoeffs_[patchi].set(facei,boundaryCoeffs_[patchi].get(facei)+delta*value); 
        }
    }
}


template<>
Foam::autoPtr<Foam::fvMatrix<Foam::scalar>::fvSolver>
Foam::fvMatrix<Foam::scalar>::solver
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solver(const dictionary& solverControls) : "
               "solver for fvMatrix<scalar>"
            << endl;
    }

    scalargpuField saveDiag(fvMatrixCache::first(level(),diag().size()),diag().size());
    saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    autoPtr<fvMatrix<scalar>::fvSolver> solverPtr
    (
        new fvMatrix<scalar>::fvSolver
        (
            *this,
            lduMatrix::solver::New
            (
                psi_.name(),
                *this,
                boundaryCoeffs_,
                internalCoeffs_,
                psi_.boundaryField().scalarInterfaces(),
                solverControls
            )
        )
    );

    diag() = saveDiag;

    return solverPtr;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::fvSolver::solve
(
    const dictionary& solverControls
)
{
    GeometricField<scalar, fvPatchField, volMesh>& psi =
        const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
        (fvMat_.psi());

    label size = fvMat_.diag().size();

    scalargpuField saveDiag(fvMatrixCache::first(fvMat_.level(),size),size);
    saveDiag = fvMat_.diag();
    fvMat_.addBoundaryDiag(fvMat_.diag(), 0);

    scalargpuField totalSource(fvMatrixCache::second(fvMat_.level(),size));
    totalSource = fvMat_.source();
    fvMat_.addBoundarySource(totalSource, false);

    // assign new solver controls
    solver_->read(solverControls);

    solverPerformance solverPerf = solver_->solve
    (
        psi.internalField(),
        totalSource
    );

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(fvMat_.mesh().comm()));
    }

    fvMat_.diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    GeometricField<scalar, fvPatchField, volMesh>& psi =
       const_cast<GeometricField<scalar, fvPatchField, volMesh>&>(psi_);

    label size = diag().size();

    scalargpuField saveDiag(fvMatrixCache::first(level(),size),size);
    saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    scalargpuField totalSource(fvMatrixCache::second(level(),size),size);
    totalSource = source_;
    addBoundarySource(totalSource, false);

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        psi.boundaryField().scalarInterfaces(),
        solverControls
    )->solve(psi.internalField(), totalSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}

namespace Foam
{

struct fvScalarMatrixResidualFunctor
{
    template<class Tuple>
    __HOST____DEVICE__ 
    scalar operator()(const scalar& source, const Tuple& t)
    {
        return source - thrust::get<0>(t)*thrust::get<1>(t);
    }
};

}

template<>
void Foam::fvMatrix<Foam::scalar>::residual(Foam::scalargpuField& tres) const
{
    scalargpuField boundaryDiag(fvMatrixCache::first(level(),psi_.size()),psi_.size());
    boundaryDiag = 0.0;
    addBoundaryDiag(boundaryDiag, 0);

    thrust::transform
    (
        source_.begin(),
        source_.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            boundaryDiag.begin(),
            psi_.internalField().begin()
        )),
        boundaryDiag.begin(),
        fvScalarMatrixResidualFunctor()
    );

    lduMatrix::residual
    (
        tres,
        psi_.internalField(),
        boundaryDiag,
        boundaryCoeffs_,
        psi_.boundaryField().scalarInterfaces(),
        0
    );

    addBoundarySource(tres);
}

template<>
Foam::tmp<Foam::scalargpuField> Foam::fvMatrix<Foam::scalar>::residual() const
{
    tmp<scalargpuField> tres(new scalargpuField(psi_.size()));

    residual(tres());

    return tres;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H() const
{
    tmp<volScalarField> tHphi
    (
        new volScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& Hphi = tHphi();

    Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().V().getField();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H1() const
{
    tmp<volScalarField> tH1
    (
        new volScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1();

    H1_.internalField() = lduMatrix::H1();
    //addBoundarySource(Hphi.internalField());

    H1_.internalField() /= psi_.mesh().V().getField();
    H1_.correctBoundaryConditions();

    return tH1;
}


// ************************************************************************* //
