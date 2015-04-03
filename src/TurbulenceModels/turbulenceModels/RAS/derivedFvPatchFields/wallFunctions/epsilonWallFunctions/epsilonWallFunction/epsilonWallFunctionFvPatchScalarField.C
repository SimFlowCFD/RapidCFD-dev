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
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar epsilonWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

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
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
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
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);
/*
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
*/
            const labelgpuList& pcells = mesh.lduAddr().patchSortCells(patchi);
            const labelgpuList& losortStart = mesh.lduAddr().patchSortStartAddr(patchi);

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
                (
                    weights.getField().begin(),
                    pcells.begin()
                ),
                addOperatorFunctor<scalar,label,scalar>()
            );
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size(), 0.0);
    epsilon_.setSize(dimensionedInternalField().size(), 0.0);

    initialised_ = true;
}


epsilonWallFunctionFvPatchScalarField&
epsilonWallFunctionFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = epsilon.boundaryField();

    const epsilonWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonWallFunctionFvPatchScalarField>(bf[patchi]);

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
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const gpuList<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

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
    const scalar* nuw;
    const scalar* nutw;
    const scalar* magGradUw;

    EpsilonCalculateGFunctor
    (
        const scalar Cmu25_,
        const scalar kappa_,
        const scalar* cornerWeights_,
        const scalar* y_,
        const scalar* k_,
        const scalar* nuw_,
        const scalar* nutw_,
        const scalar* magGradUw_
    ):
        Cmu25(Cmu25_),
        kappa(kappa_),
        cornerWeights(cornerWeights_),
        y(y_),
        k(k_),
        nuw(nuw_),
        nutw(nutw_),
        magGradUw(magGradUw_)
    {}

    __HOST____DEVICE__
    scalar operator()(const label& cellI, const label& faceI)
    {
        scalar w = cornerWeights[faceI];

        return
            w
           *(nutw[faceI] + nuw[faceI])
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
    const label patchi = patch.index();

    const scalargpuField& y = turbulence.y()[patchi];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<scalargpuField> tnuw = turbulence.nu(patchi);
    const scalargpuField& nuw = tnuw();

    const tmp<scalargpuField> tnutw = turbulence.nut(patchi);
    const scalargpuField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchi];

    const scalargpuField magGradUw(mag(Uw.snGrad()));

    matrixPatchOperation
    (
        patchi,
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
        patchi,
        G,
        patch.boundaryMesh().mesh().lduAddr(),
        EpsilonCalculateGFunctor
        (
            Cmu25,
            kappa_,
            cornerWeights.data(),
            y.data(),
            k.getField().data(),
            nuw.data(),
            nutw.data(),
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

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());
/*
    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        epsilon[celli] = epsilon0[celli];
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

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalargpuField& G0 = this->G();
    const scalargpuField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& epsilon = const_cast<FieldType&>(dimensionedInternalField());

    scalargpuField& epsilonf = *this;

    // only set the values if the weights are > tolerance
/*
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
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
            (epsilon.getField().begin(),patch().faceCells().begin()
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

    labelgpuList constraintCells(weights.size());
    scalargpuList constraintEpsilon(weights.size());
    const labelgpuList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon
        = dimensionedInternalField();

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
    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintEpsilon.append(epsilon[celli]);
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

} // End namespace Foam

// ************************************************************************* //
