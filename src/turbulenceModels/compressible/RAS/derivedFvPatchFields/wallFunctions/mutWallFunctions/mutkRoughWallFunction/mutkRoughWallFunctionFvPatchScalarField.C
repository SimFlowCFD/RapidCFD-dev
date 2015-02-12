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

#include "mutkRoughWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

__HOST____DEVICE__
scalar fnRough
(
    const scalar KsPlus,
    const scalar Cs
)
{
    // Return fn based on non-dimensional roughness height

    if (KsPlus < 90.0)
    {
        return pow
        (
            (KsPlus - 2.25)/87.75 + Cs*KsPlus,
            sin(0.4258*(log(KsPlus) - 0.811))
        );
    }
    else
    {
        return (1.0 + Cs*KsPlus);
    }
}

struct mutkRoughCalcMutFunctor
{
    const scalar Cmu25;
    const scalar kappa;
    const scalar E;

    mutkRoughCalcMutFunctor(const scalar Cmu25_,const scalar kappa_,const scalar E_)
                           :Cmu25(Cmu25_),kappa(kappa_),E(E_){}
    __HOST____DEVICE__
    scalar operator () (const scalar& y, const thrust::tuple<scalar,scalar,scalar,scalar,scalar,scalar>& t)
    {
        scalar k = thrust::get<0>(t);
        scalar muw = thrust::get<1>(t);
        scalar mutw = thrust::get<2>(t);
        scalar rhow = thrust::get<3>(t);
        scalar Ks = thrust::get<4>(t);
        scalar Cs = thrust::get<5>(t);
		
        scalar uStar = Cmu25*sqrt(k);
        scalar yPlus = uStar*y/(muw/rhow);
        scalar KsPlus = uStar*Ks/(muw/rhow);

        scalar Edash = E;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs);
        }

        scalar limitingMutw = max(mutw, muw);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        return
            max
            (
                min
                (
                    muw
                   *(yPlus*kappa/log(max(Edash*yPlus, 1+1e-4)) - 1),
                    2*limitingMutw
                ), 0.5*limitingMutw
            );
    }
};


tmp<scalargpuField> mutkRoughWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalargpuField& y = turbModel.y()[patchi];
    const scalargpuField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalargpuField& muw = turbModel.mu().boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalargpuField> tmutw(new scalargpuField(*this));
    scalargpuField& mutw = tmutw();
    
    thrust::transform
    (
        y.begin(),
        y.end(),
        thrust::make_zip_iterator(thrust::make_tuple
        (
            thrust::make_permutation_iterator
            (
                k.getField().begin(),
                patch().faceCells().begin()
            ),
            muw.begin(),
            mutw.begin(),
            rhow.begin(),
            Ks_.begin(),
            Cs_.begin()
        )),
        mutw.begin(),
        mutkRoughCalcMutFunctor(Cmu25,kappa_,E_)
    );

/*
    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uStar = Cmu25*sqrt(k[faceCellI]);
        scalar yPlus = uStar*y[faceI]/(muw[faceI]/rhow[faceI]);
        scalar KsPlus = uStar*Ks_[faceI]/(muw[faceI]/rhow[faceI]);

        scalar Edash = E_;
        if (KsPlus > 2.25)
        {
            Edash /= fnRough(KsPlus, Cs_[faceI]);
        }

        scalar limitingMutw = max(mutw[faceI], muw[faceI]);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        mutw[faceI] =
            max
            (
                min
                (
                    muw[faceI]
                   *(yPlus*kappa_/log(max(Edash*yPlus, 1+1e-4)) - 1),
                    2*limitingMutw
                ), 0.5*limitingMutw
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", mutw = " << mutw[faceI]
                << endl;
        }
    }
*/
    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), 0.0),
    Cs_(p.size(), 0.0)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size())
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


mutkRoughWallFunctionFvPatchScalarField::mutkRoughWallFunctionFvPatchScalarField
(
    const mutkRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mutkWallFunctionFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void mutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelgpuList& addr
)
{
    mutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const mutkRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const mutkRoughWallFunctionFvPatchScalarField>(ptf);

    Cs_.rmap(nrwfpsf.Cs_, addr);
    Ks_.rmap(nrwfpsf.Ks_, addr);
}


void mutkRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutkRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
