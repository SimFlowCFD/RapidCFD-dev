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

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar epsilonWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

struct epsilonWallFunctionGraterThanToleranceFunctor : public std::unary_function<scalar,bool>
{
    const scalar tolerance_;

    epsilonWallFunctionGraterThanToleranceFunctor
    (
        scalar tolerance
    ): 
        tolerance_(tolerance)
    {}

    __HOST____DEVICE__
    bool operator () (const scalar& s)
    {
        return s > tolerance_;
    }
};

void epsilonWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("epsilonWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void epsilonWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


void epsilonWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchI)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            if (master == -1)
            {
                master = patchI;
            }

            epf.master() = master;
        }
    }
}

struct epsilonWallFunctionFvPatchScalarFieldCreateWeightsFunctor: public std::unary_function<thrust::tuple<label,label>,label>
{
    __HOST____DEVICE__
    label operator()(const thrust::tuple<label,label>& t)
    {
        return thrust::get<1>(t) - thrust::get<0>(t);
    }
};

void epsilonWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

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

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchI)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchI]))
        {
/*            epsilonPatches.append(patchI);

            const labelUList& faceCells = bf[patchI].patch().faceCells();
            forAll(faceCells, i)
            {
                label cellI = faceCells[i];
                weights[cellI]++;
            }
*/            
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
                    epsilonWallFunctionFvPatchScalarFieldCreateWeightsFunctor()
                ),
                thrust::make_permutation_iterator
                (weights.getField().begin(),pcells.begin()
                ),
                addOperatorFunctor<scalar,label,scalar>()
            );
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchI = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchI];
        cornerWeights_[patchI] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size());
    epsilon_.setSize(dimensionedInternalField().size());

    initialised_ = true;
}


epsilonWallFunctionFvPatchScalarField&
epsilonWallFunctionFvPatchScalarField::epsilonPatch(const label patchI)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    const epsilonWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonWallFunctionFvPatchScalarField>(bf[patchI]);

    return const_cast<epsilonWallFunctionFvPatchScalarField&>(epf);
}


void epsilonWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalargpuField& G0,
    scalargpuField& epsilon0
)
{
    // accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            const gpuList<scalar>& w = cornerWeights_[patchI];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchI);

            epf == scalargpuField(epsilon0, epf.patch().faceCells());
        }
    }
}

    
struct EpsilonCalculateEpsilonFunctor : public std::unary_function<label,scalar>
{
    const scalar Cmu75;
    const scalar kappa;
    const scalar* cornerWeights;
    const scalar* y;
    const scalar* k;

    EpsilonCalculateEpsilonFunctor
    (
        const scalar Cmu75_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_
    ):

        Cmu75(Cmu75_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        return w*Cmu75*pow(k[cellI], 1.5)/(kappa*y[faceI]);
    }
};

struct EpsilonCalculateGFunctor : public std::unary_function<label,scalar>
{
    const scalar Cmu25;
    const scalar kappa;
    const scalar* cornerWeights;
    const scalar* y;
    const scalar* k;
    const scalar* muw;
    const scalar* mutw;
    const scalar* magGradUw;

    EpsilonCalculateGFunctor
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


void epsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const gpuList<scalar>& cornerWeights,
    const fvPatch& patch,
    scalargpuField& G,
    scalargpuField& epsilon
)
{
    const label patchI = patch.index();

    const scalargpuField& y = turbulence.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalargpuField& muw = tmu().boundaryField()[patchI];

    const tmp<volScalarField> tmut = turbulence.mut();
    const volScalarField& mut = tmut();
    const scalargpuField& mutw = mut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalargpuField magGradUw(mag(Uw.snGrad()));

    // Set epsilon and G
                                               
    matrixPatchOperation
    (
        patchI,
        epsilon,
        patch.boundaryMesh().mesh().lduAddr(),
        EpsilonCalculateEpsilonFunctor
        (
            Cmu75,
            kappa_,
            cornerWeights.data(),
            y.data(),
            k.getField().data()
        )
    );
						                        
	
    matrixPatchOperation
    (
        patchI,
        G,
        patch.boundaryMesh().mesh().lduAddr(),
        EpsilonCalculateGFunctor
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

epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
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
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    G_(),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalargpuField& epsilonWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


scalargpuField& epsilonWallFunctionFvPatchScalarField::epsilon(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void epsilonWallFunctionFvPatchScalarField::updateCoeffs()
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
        calculateTurbulenceFields(turbulence, G(true), epsilon(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());
/*
    forAll(*this, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        G[cellI] = G0[cellI];
        epsilon[cellI] = epsilon0[cellI];
    }
*/
    thrust::copy(thrust::make_permutation_iterator(G0.begin(),patch().faceCells().begin()),
                 thrust::make_permutation_iterator(G0.begin(),patch().faceCells().end()),
                 thrust::make_permutation_iterator(G.getField().begin(),patch().faceCells().begin()));
                 
    thrust::copy(thrust::make_permutation_iterator(epsilon0.begin(),patch().faceCells().begin()),
                 thrust::make_permutation_iterator(epsilon0.begin(),patch().faceCells().end()),
                 thrust::make_permutation_iterator(epsilon.getField().begin(),patch().faceCells().begin()));

    fvPatchField<scalar>::updateCoeffs();
}

struct epsilonWallFunctionFvPatchScalarFieldupdateCoeffsFunctor
{
    __HOST____DEVICE__
    scalar operator () (const scalar& w,const thrust::tuple<scalar,scalar>& t)
    {
        scalar G = thrust::get<0>(t);
        scalar G0 = thrust::get<1>(t);
		
        return (1.0-w)*G+w*G0;
    }	
};

void epsilonWallFunctionFvPatchScalarField::updateCoeffs
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
        calculateTurbulenceFields(turbulence, G(true), epsilon(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());

    scalargpuField& epsilonf = *this;

    // only set the values if the weights are > tolerance
    /*
    forAll(weights, faceI)
    {
        scalar w = weights[faceI];

        if (w > tolerance_)
        {
            label cellI = patch().faceCells()[faceI];

            G[cellI] = (1.0 - w)*G[cellI] + w*G0[cellI];
            epsilon[cellI] = (1.0 - w)*epsilon[cellI] + w*epsilon0[cellI];
            epsilonf[faceI] = epsilon[cellI];
        }
    }
    */

    thrust::transform
    (
        weights.begin(),
        weights.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator
            (
                G.getField().begin(),
                patch().faceCells().begin()
            ),
            thrust::make_permutation_iterator
            (
                G0.begin(),
                patch().faceCells().begin()
            )
        )),
        thrust::make_permutation_iterator
        (
            G.getField().begin(),
            patch().faceCells().begin()
        ),
        epsilonWallFunctionFvPatchScalarFieldupdateCoeffsFunctor()
    );
                     
    thrust::transform
    (
        weights.begin(),
        weights.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator
            (
                epsilon.getField().begin(),
                patch().faceCells().begin()
            ),
            thrust::make_permutation_iterator
            (
                epsilon0.begin(),
                patch().faceCells().begin()
            )
        )),
        thrust::make_permutation_iterator
        (
            epsilon.getField().begin(),
            patch().faceCells().begin()
        ),
        epsilonWallFunctionFvPatchScalarFieldupdateCoeffsFunctor()
    );
                     
                     
    thrust::copy
    (
        thrust::make_permutation_iterator
        (
            epsilon.getField().begin(),
            patch().faceCells().begin()
        ),
        thrust::make_permutation_iterator
        (
            epsilon.getField().begin(),
            patch().faceCells().end()
        ),
        epsilonf.begin()
    );

    fvPatchField<scalar>::updateCoeffs();
}


void epsilonWallFunctionFvPatchScalarField::manipulateMatrix
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


void epsilonWallFunctionFvPatchScalarField::manipulateMatrix
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
    scalargpuList constraintEpsilon(weights.size());
    const labelgpuList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon =
        dimensionedInternalField();

    typename gpuList<label>::iterator end = 
        thrust::copy_if
        (
            faceCells.begin(),
            faceCells.end(),
            weights.begin(),
            constraintCells.begin(),
            epsilonWallFunctionGraterThanToleranceFunctor(tolerance_)
        );
                               
    label nConstrainedCells = end - constraintCells.begin();
    
    constraintCells.setSize(nConstrainedCells);
    
    thrust::copy_if
    (
        thrust::make_permutation_iterator
        (
            epsilon.getField().begin(),
            faceCells.begin()
        ),
        thrust::make_permutation_iterator
        (
            epsilon.getField().begin(),
            faceCells.end()
        ),
        weights.begin(),
        constraintEpsilon.begin(),
        epsilonWallFunctionGraterThanToleranceFunctor(tolerance_)
    );
                               
    constraintEpsilon.setSize(nConstrainedCells);
/*
    forAll(weights, faceI)
    {
        // only set the values if the weights are > tolerance
        if (weights[faceI] > tolerance_)
        {
            nConstrainedCells++;

            label cellI = faceCells[faceI];

            constraintCells.append(cellI);
            constraintEpsilon.append(epsilon[cellI]);
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
                        constraintEpsilon.begin());
                        
    matrix.setValues
    (
        constraintCells,
        scalargpuField(constraintEpsilon.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void epsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
