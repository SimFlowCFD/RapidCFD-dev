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

#include "omegaWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "mutWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

inline scalar yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}

struct omegaWallFunctionGraterThanToleranceFunctor : public std::unary_function<scalar,bool>
{
    const scalar tolerance_;
    omegaWallFunctionGraterThanToleranceFunctor(scalar tolerance): tolerance_(tolerance){}
    bool operator () (const scalar& s)
    {
        return s > tolerance_;
    }
};

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("omegaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void omegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}


void omegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchI)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            if (master == -1)
            {
                master = patchI;
            }

            opf.master() = master;
        }
    }
}

struct omegaWallFunctionFvPatchScalarFieldCreateWeightsFunctor: public std::unary_function<thrust::tuple<label,label>,label>
{
    __HOST____DEVICE__
    label operator()(const thrust::tuple<label,label>& t)
    {
        return thrust::get<1>(t) - thrust::get<0>(t);
    }
};

void omegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchI)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            omegaPatches.append(patchI);

            const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchI);
            const labelgpuList& losortStart = mesh.lduAddr().patchSortStartAddr(patchI);
            
            thrust::transform
            (
                thrust::make_permutation_iterator
                (
                    weights.getField().begin(),
                    pcells.begin()
                ),
                thrust::make_permutation_iterator
                (
                    weights.getField().begin(),
                    pcells.end()
                ),
                thrust::make_transform_iterator
                (
                    thrust::make_zip_iterator(thrust::make_tuple
                    (
                        losortStart.begin(),
                        losortStart.begin()+1
                    )),
                    omegaWallFunctionFvPatchScalarFieldCreateWeightsFunctor()
                ),
                thrust::make_permutation_iterator(weights.getField().begin(),pcells.begin()),
                addOperatorFunctor<scalar,label,scalar>()
            );           
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchI = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchI];
        cornerWeights_[patchI] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size());
    omega_.setSize(dimensionedInternalField().size());

    initialised_ = true;
}


omegaWallFunctionFvPatchScalarField&
omegaWallFunctionFvPatchScalarField::omegaPatch(const label patchI)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const omegaWallFunctionFvPatchScalarField& opf =
        refCast<const omegaWallFunctionFvPatchScalarField>(bf[patchI]);

    return const_cast<omegaWallFunctionFvPatchScalarField&>(opf);
}


void omegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalargpuField& G0,
    scalargpuField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            const gpuList<scalar>& w = cornerWeights_[patchI];

            opf.calculate(turbulence, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            opf == scalargpuField(omega0, opf.patch().faceCells());
        }
    }
}

    
struct OmegaCalculateOmegaFunctor : public std::unary_function<label,scalar>
{
    const scalar Cmu25;
    const scalar kappa;
    const scalar beta1;
    const scalar* cornerWeights;
    const scalar* y;
    const scalar* k;
    const scalar* muw;
    const scalar* rhow;

    OmegaCalculateOmegaFunctor
    (
        const scalar Cmu25_,
        const scalar kappa_,
        const scalar beta1_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* muw_,
        const scalar* rhow_
    ):
        Cmu25(Cmu25_),
        kappa(kappa_),
        beta1(beta1_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        muw(muw_),
        rhow(rhow_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI,const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        scalar omegaVis = 6.0*muw[faceI]/(rhow[faceI]*beta1*sqr(y[faceI]));

        scalar omegaLog = sqrt(k[cellI])/(Cmu25*kappa*y[faceI]);

        return w*sqrt(sqr(omegaVis) + sqr(omegaLog));
    }

};

struct OmegaCalculateGFunctor : public std::unary_function<label,scalar>
{
    const scalar Cmu25;
    const scalar kappa;
    const scalar* cornerWeights;
    const scalar* y;
    const scalar* k;
    const scalar* muw;
    const scalar* mutw;
    const scalar* magGradUw;

    OmegaCalculateGFunctor
    (
        const scalar Cmu25_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* muw_,
        const scalar* mutw_,
        const scalar* magGradUw_
    ):
        Cmu25(Cmu25_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        muw(muw_),
        mutw(mutw_),
        magGradUw(magGradUw_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        return
            w
           *(mutw[faceI] + muw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[cellI])
           /(kappa*y[faceI]);
    }
};

void omegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const gpuList<scalar>& cornerWeights,
    const fvPatch& patch,
    scalargpuField& G,
    scalargpuField& omega
)
{
    const label patchI = patch.index();

    const scalargpuField& y = turbulence.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const scalargpuField& rhow = turbulence.rho().boundaryField()[patchI];

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalargpuField& muw = tmu().boundaryField()[patchI];

    const tmp<volScalarField> tmut = turbulence.mut();
    const volScalarField& mut = tmut();
    const scalargpuField& mutw = mut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalargpuField magGradUw(mag(Uw.snGrad()));

    // Set omega and G	
                                               
    matrixPatchOperation
    (
        patchI,
        omega,
        patch.boundaryMesh().mesh().lduAddr(),
        OmegaCalculateOmegaFunctor
        (
            Cmu25,
            kappa_,
            beta1_,
            cornerWeights.data(),
            y.data(),
            k.getField().data(),
            muw.data(),
            rhow.data()
        )
    );
						                        	                        
	
    matrixPatchOperation
    (
        patchI,
        G,
        patch.boundaryMesh().mesh().lduAddr(),
        OmegaCalculateGFunctor
        (
            Cmu25,
            kappa_,
            cornerWeights.data(),
            y.data(),
            k.getField().data(),
            muw.data(),
            mutw.data(),
            magGradUw.data()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075),
    yPlusLam_(yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    yPlusLam_(yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalargpuField& omegaWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalargpuField& omegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());
/*
    forAll(*this, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        G[cellI] = G0[cellI];
        omega[cellI] = omega0[cellI];
    }
*/
    thrust::copy(thrust::make_permutation_iterator(G0.begin(),patch().faceCells().begin()),
                 thrust::make_permutation_iterator(G0.begin(),patch().faceCells().end()),
                 thrust::make_permutation_iterator(G.getField().begin(),patch().faceCells().begin()));
                 
    thrust::copy(thrust::make_permutation_iterator(omega0.begin(),patch().faceCells().begin()),
                 thrust::make_permutation_iterator(omega0.begin(),patch().faceCells().end()),
                 thrust::make_permutation_iterator(omega.getField().begin(),patch().faceCells().begin()));

    fvPatchField<scalar>::updateCoeffs();
}

struct omegaWallFunctionFvPatchScalarFieldupdateCoeffsFunctor
{
    __HOST____DEVICE__
    scalar operator () (const scalar& w,const thrust::tuple<scalar,scalar>& t)
    {
        scalar G = thrust::get<0>(t);
        scalar G0 = thrust::get<1>(t);
		
        return (1.0-w)*G+w*G0;
    }	
	
};

void omegaWallFunctionFvPatchScalarField::updateCoeffs
(
    const scalargpuField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());

    scalargpuField& omegaf = *this;

    // only set the values if the weights are > tolerance_
    /*
    forAll(weights, faceI)
    {
        scalar w = weights[faceI];

        if (w > tolerance_)
        {
            label cellI = patch().faceCells()[faceI];

            G[cellI] = (1.0 - w)*G[cellI] + w*G0[cellI];
            omega[cellI] = (1.0 - w)*omega[cellI] + w*omega0[cellI];
            omegaf[faceI] = omega[cellI];
        }
    }
    */

    thrust::transform
    (
        weights.begin(),
        weights.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
             thrust::make_permutation_iterator(G.getField().begin(),  patch().faceCells().begin()),
             thrust::make_permutation_iterator(G0.begin(),            patch().faceCells().begin())
        )),
        thrust::make_permutation_iterator(G.getField().begin(),       patch().faceCells().begin()),
        omegaWallFunctionFvPatchScalarFieldupdateCoeffsFunctor()
    );
                     
    thrust::transform
    (
        weights.begin(),
        weights.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
             thrust::make_permutation_iterator(omega.getField().begin(),   patch().faceCells().begin()),
             thrust::make_permutation_iterator(omega0.begin(),             patch().faceCells().begin())
        )),
        thrust::make_permutation_iterator(omega.getField().begin(),patch().faceCells().begin()),
        omegaWallFunctionFvPatchScalarFieldupdateCoeffsFunctor()
    );
                     
                     
    thrust::copy
    (
        thrust::make_permutation_iterator(omega.getField().begin(),patch().faceCells().begin()),
        thrust::make_permutation_iterator(omega.getField().begin(),patch().faceCells().end()),
        omegaf.begin()
    );

    fvPatchField<scalar>::updateCoeffs();
}


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const gpuField<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    // filter weights so that we only apply the constraint where the
    // weight > SMALL
    labelgpuList constraintCells(weights.size());
    scalargpuList constraintOmega(weights.size());
    const labelgpuList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega =
        dimensionedInternalField();

    typename gpuList<label>::iterator end = 
        thrust::copy_if
        (
            faceCells.begin(),
            faceCells.end(),
            weights.begin(),
            constraintCells.begin(),
            omegaWallFunctionGraterThanToleranceFunctor(tolerance_)
        );
                               
    label nConstrainedCells = end - constraintCells.begin();
    
    constraintCells.setSize(nConstrainedCells);
    
    thrust::copy_if
    (
        thrust::make_permutation_iterator(omega.getField().begin(),faceCells.begin()),
        thrust::make_permutation_iterator(omega.getField().begin(),faceCells.end()),
        weights.begin(),
        constraintOmega.begin(),
        omegaWallFunctionGraterThanToleranceFunctor(tolerance_)
    );
                               
    constraintOmega.setSize(nConstrainedCells);
/*
    forAll(weights, faceI)
    {
        // only set the values if the weights are > tolerance
        if (weights[faceI] > tolerance_)
        {
            nConstrainedCells++;

            label cellI = faceCells[faceI];

            constraintCells.append(cellI);
            constraintomega.append(omega[cellI]);
        }
    }
*/
    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }
    
    thrust::sort_by_key(constraintCells.begin(),
                        constraintCells.end(),
                        constraintOmega.begin());

    matrix.setValues
    (
        constraintCells,
        scalargpuField(constraintOmega.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
