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

#include "LduMatrix.H"
#include "diagTensorField.H"
#include "fvMatrixCache.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
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


template<class Type>
Foam::solverPerformance Foam::fvMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solve(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    label maxIter = -1;
    if (solverControls.readIfPresent("maxIter", maxIter))
    {
        if (maxIter == 0)
        {
            return solverPerformance();
        }
    }

    word type(solverControls.lookupOrDefault<word>("type", "segregated"));

    if (type == "segregated")
    {
        return solveSegregated(solverControls);
    }
    else if (type == "coupled")
    {
        return solveCoupled(solverControls);
    }
    else
    {
        FatalIOErrorIn
        (
            "fvMatrix<Type>::solve(const dictionary& solverControls)",
            solverControls
        )   << "Unknown type " << type
            << "; currently supported solver types are segregated and coupled"
            << exit(FatalIOError);

        return solverPerformance();
    }
}


template<class Type>
Foam::solverPerformance Foam::fvMatrix<Type>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    GeometricField<Type, fvPatchField, volMesh>& psi =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    solverPerformance solverPerfVec
    (
        "fvMatrix<Type>::solveSegregated",
        psi.name()
    );

   label size = diag().size();

    scalargpuField saveDiag(fvMatrixCache::first(level(),size),size);
    saveDiag = diag();

    gpuField<Type> source(source_);

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implict part by updateMatrixInterfaces within
    // the component loop.
    addBoundarySource(source);

    typename Type::labelType validComponents
    (
        pow
        (
            psi.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // copy field and source

        scalargpuField psiCmpt(fvMatrixCache::second(level(),size),size);
        component(psiCmpt,psi.internalField(),cmpt);
        addBoundaryDiag(diag(), cmpt);

        scalargpuField sourceCmpt(fvMatrixCache::third(level(),size),size);
        component(sourceCmpt,source,cmpt);

        FieldField<gpuField, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<gpuField, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        lduInterfaceFieldPtrsList interfaces =
            psi.boundaryField().scalarInterfaces();

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        initMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        if (solverPerformance::debug)
        {
            solverPerf.print(Info.masterStream(this->mesh().comm()));
        }

        solverPerfVec = max(solverPerfVec, solverPerf);
        solverPerfVec.solverName() = solverPerf.solverName();

        psi.internalField().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerfVec);

    return solverPerfVec;
}


template<class Type>
Foam::solverPerformance Foam::fvMatrix<Type>::solveCoupled
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solveCoupled"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    GeometricField<Type, fvPatchField, volMesh>& psi =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    LduMatrix<Type, scalar, scalar> coupledMatrix(psi.mesh());
    coupledMatrix.diag() = diag();
    coupledMatrix.upper() = upper();
    coupledMatrix.lower() = lower();
    coupledMatrix.source() = source();

    addBoundaryDiag(coupledMatrix.diag(), 0);
    addBoundarySource(coupledMatrix.source(), false);

    coupledMatrix.interfaces() = psi.boundaryField().interfaces();
    coupledMatrix.interfacesUpper() = boundaryCoeffs().component(0);
    coupledMatrix.interfacesLower() = internalCoeffs().component(0);

    autoPtr<typename LduMatrix<Type, scalar, scalar>::solver>
    coupledMatrixSolver
    (
        LduMatrix<Type, scalar, scalar>::solver::New
        (
            psi.name(),
            coupledMatrix,
            solverControls
        )
    );

    SolverPerformance<Type> solverPerf
    (
        coupledMatrixSolver->solve(psi.getField())
    );

    if (SolverPerformance<Type>::debug)
    {
        solverPerf.print(Info);
    }

    psi.correctBoundaryConditions();

    // psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerformance();
}


template<class Type>
Foam::autoPtr<typename Foam::fvMatrix<Type>::fvSolver>
Foam::fvMatrix<Type>::solver()
{
    return solver
    (
        psi_.mesh().solverDict
        (
            psi_.select
            (
                psi_.mesh().data::template lookupOrDefault<bool>
                ("finalIteration", false)
            )
        )
    );
}


template<class Type>
Foam::solverPerformance Foam::fvMatrix<Type>::fvSolver::solve()
{
    return solve
    (
        fvMat_.psi_.mesh().solverDict
        (
            fvMat_.psi_.select
            (
                fvMat_.psi_.mesh().data::template lookupOrDefault<bool>
                ("finalIteration", false)
            )
        )
    );
}


template<class Type>
Foam::solverPerformance Foam::fvMatrix<Type>::solve()
{
    return solve
    (
        psi_.mesh().solverDict
        (
            psi_.select
            (
                psi_.mesh().data::template lookupOrDefault<bool>
                ("finalIteration", false)
            )
        )
    );
}

template<class Type>
void Foam::fvMatrix<Type>::residual(Foam::gpuField<Type>& res) const
{
    res = source_;

    addBoundarySource(res);

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        label pSize = psi_.size();

        scalargpuField psiCmpt(fvMatrixCache::first(level(),pSize),pSize);
        component(psiCmpt,psi_.internalField(),cmpt);

        scalargpuField boundaryDiagCmpt(fvMatrixCache::second(level(),pSize),pSize);
        boundaryDiagCmpt = 0.0;

        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<gpuField, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                psi_.boundaryField().scalarInterfaces(),
                cmpt
            )
        );
    }

}


template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::fvMatrix<Type>::residual() const
{
    tmp<gpuField<Type> > tres(new gpuField<Type>(source_));
    gpuField<Type>& res = tres();

    residual(res);

    return tres;
}


// ************************************************************************* //
