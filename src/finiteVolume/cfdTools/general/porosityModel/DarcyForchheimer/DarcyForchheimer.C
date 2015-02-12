/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "DarcyForchheimer.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::DarcyForchheimer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    dXYZ_(coeffs_.lookup("d")),
    fXYZ_(coeffs_.lookup("f")),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu"))
{
    adjustNegativeResistance(dXYZ_);
    adjustNegativeResistance(fXYZ_);

    calcTranformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::DarcyForchheimer::calcTranformModelData()
{
    if (coordSys_.R().uniform())
    {
        forAll (cellZoneIDs_, zoneI)
        {
            D_[zoneI].setSize(1);
            F_[zoneI].setSize(1);
            
            tensor dr = tensor::zero;
            dr.xx() = dXYZ_.value().x();
            dr.yy() = dXYZ_.value().y();
            dr.zz() = dXYZ_.value().z();
            
            D_[zoneI].set(0,coordSys_.R().transformTensor(dr));

            // leading 0.5 is from 1/2*rho
            tensor fr = tensor::zero;
            fr.xx() = 0.5*fXYZ_.value().x();
            fr.yy() = 0.5*fXYZ_.value().y();
            fr.zz() = 0.5*fXYZ_.value().z();
            F_[zoneI].set(0,coordSys_.R().transformTensor(fr));
        }
    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();

            D_[zoneI].setSize(cells.size());
            F_[zoneI].setSize(cells.size());

            tensor dr = tensor::zero;
            dr.xx() = dXYZ_.value().x();
            dr.yy() = dXYZ_.value().y();
            dr.zz() = dXYZ_.value().z();
            
            D_[zoneI].operator=(dr);

            // leading 0.5 is from 1/2*rho
            tensor fr = tensor::zero;
            fr.xx() = 0.5*fXYZ_.value().x();
            fr.yy() = 0.5*fXYZ_.value().y();
            fr.zz() = 0.5*fXYZ_.value().z();
            F_[zoneI].operator=(fr);

            const coordinateRotation& R = coordSys_.R(mesh_, cells);

            D_[zoneI] = R.transformTensor(D_[zoneI], cells);
            F_[zoneI] = R.transformTensor(F_[zoneI], cells);
            

        }
    }
/*
    if (debug && mesh_.time().outputTime())
    {
        volTensorField Dout
        (
            IOobject
            (
                typeName + ":D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", dXYZ_.dimensions(), tensor::zero)
        );
        volTensorField Fout
        (
            IOobject
            (
                typeName + ":F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", fXYZ_.dimensions(), tensor::zero)
        );

        UIndirectList<tensor>(Dout, mesh_.cellZones()[cellZoneIDs_[0]]) = D_[0];
        UIndirectList<tensor>(Fout, mesh_.cellZones()[cellZoneIDs_[0]]) = F_[0];

        Dout.write();
        Fout.write();
    }
*/
}


void Foam::porosityModels::DarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorgpuField& force
) const
{
    scalargpuField Udiag(U.size(), 0.0);
    vectorgpuField Usource(U.size(), vector::zero);
    const scalargpuField& V = mesh_.V().getField();

    apply(Udiag, Usource, V, rho.getField(), mu.getField(), U.getField());

    force = Udiag*U.getField() - Usource;
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalargpuField& V = mesh_.V().getField();
    scalargpuField& Udiag = UEqn.diag();
    vectorgpuField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho.getField(), mu.getField(), U.getField());
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho.getField(), (rho*nu)().getField(), U.getField());
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu.getField(), U.getField());
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), (mu/rho)().getField(), U.getField());
        }
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorgpuField& U = UEqn.psi().getField();
    const scalargpuField& V = mesh_.V().getField();
    scalargpuField& Udiag = UEqn.diag();
    vectorgpuField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho.getField(), mu.getField(), U);
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU.getField(), rho.getField(), mu.getField(), U.getField());
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU.getField(), geometricOneField(), nu.getField(), U.getField());
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU.getField(), geometricOneField(), (mu/rho)().getField(), U.getField());
        }
    }
}


bool Foam::porosityModels::DarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
