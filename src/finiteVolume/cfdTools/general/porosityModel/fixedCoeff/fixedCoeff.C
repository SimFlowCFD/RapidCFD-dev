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
#include "fixedCoeff.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porosityModel, fixedCoeff, mesh);
    }
    
    __HOST____DEVICE__
    inline tensor calcCd(const tensor& alpha, const tensor& beta, const scalar& rho, const vector& U)
    {
        return rho*(alpha + beta*mag(U));
    }
	
    
    struct fixeCoeffAUFunctor
    {
        const scalar rho;
        fixeCoeffAUFunctor(scalar _rho): rho(_rho) {}
        __HOST____DEVICE__
        tensor operator () (const tensor& AU, const thrust::tuple<tensor,tensor,vector>& t)
        {
            return AU + calcCd(thrust::get<0>(t),thrust::get<1>(t),rho,thrust::get<2>(t));
        }
    };
	
    struct fixeCoeffConstAUFunctor
    {
        const tensor alpha;
        const tensor beta;
        const scalar rho;
        fixeCoeffConstAUFunctor(tensor _alpha, tensor _beta, scalar _rho): alpha(_alpha), beta(_beta), rho(_rho) {}
        __HOST____DEVICE__
        tensor operator () (const tensor& AU, const thrust::tuple<vector>& t)
        {
            return AU + calcCd(alpha,beta,rho,thrust::get<0>(t));
        }
     };
	
	
    __HOST____DEVICE__
    inline scalar calcUdiag(const scalar& V,const scalar& isoCd)
    {
        return V*isoCd;
    }

    __HOST____DEVICE__
    inline vector calcUsource(const scalar& V,const scalar& isoCd, const tensor& Cd,const vector& U)
    {
        return V*((Cd - tensor(1,0,0,0,1,0,0,0,1)*isoCd) & U);
    }

    __HOST____DEVICE__
    inline thrust::tuple<scalar,vector> DFapply(const tensor& alpha, const tensor& beta,const scalar& rho,const vector& U,const scalar& V)
    {
        const tensor Cd = calcCd(alpha,beta,rho,U);
        const scalar isoCd = tr(Cd);
		
        scalar Udiag = calcUdiag(V,isoCd);
        vector Usource = calcUsource(V,isoCd,Cd,U);
		
        return thrust::make_tuple(Udiag,Usource);
    }

    struct fixedCoeffFunctor
    {
        const scalar rho;
        fixedCoeffFunctor(scalar _rho): rho(_rho) {}
        __HOST____DEVICE__
        thrust::tuple<scalar,vector> operator ()(const thrust::tuple<scalar,vector>& in,const thrust::tuple<tensor,tensor,vector,scalar>& t)
        {
            thrust::tuple<scalar,vector>  out =  
                DFapply
                (
                    thrust::get<0>(t),
                    thrust::get<1>(t),
                    rho,
                    thrust::get<2>(t),
                    thrust::get<3>(t)
                );
                    

            return thrust::make_tuple
            (
                thrust::get<0>(in)+thrust::get<0>(out),
                thrust::get<1>(in)+thrust::get<1>(out)
            );

        }
    };

    struct fixedCoeffConstFunctor
    {
        const tensor alpha;
        const tensor beta;
        const scalar rho;
        fixedCoeffConstFunctor(tensor alpha_, tensor beta_, scalar _rho): alpha(alpha_), beta(beta_), rho(_rho) {}
        __HOST____DEVICE__
        thrust::tuple<scalar,vector> operator ()(const thrust::tuple<scalar,vector>& in,const thrust::tuple<vector,scalar>& t)
        {
            thrust::tuple<scalar,vector>  out =  
                DFapply
                (
                    alpha,
                    beta,
                    rho,
                    thrust::get<0>(t),
                    thrust::get<1>(t)
                );

            return thrust::make_tuple
                (
                    thrust::get<0>(in)+thrust::get<0>(out),
                    thrust::get<1>(in)+thrust::get<1>(out)
                );
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::apply
(
    scalargpuField& Udiag,
    vectorgpuField& Usource,
    const scalargpuField& V,
    const vectorgpuField& U,
    const scalar rho
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorgpuField& alphaZones = alpha_[zoneI];
        const tensorgpuField& betaZones = beta_[zoneI];

        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();

        if(alphaZones.size() == 1)
        {
            tensor alpha = alphaZones.get(0);
            tensor beta = betaZones.get(0);

            thrust::transform
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.end()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.end()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        U.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        V.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.begin()
                    )
                )),
                fixedCoeffConstFunctor(alpha,beta,rho)
            );
        }
        else
        {
            thrust::transform
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.end()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.end()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    alphaZones.begin(),
                    betaZones.begin(),
                    thrust::make_permutation_iterator
                    (
                        U.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        V.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        Udiag.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        Usource.begin(),
                        cells.begin()
                    )
                )),
                fixedCoeffFunctor(rho)
            );
        }
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorgpuField& AU,
    const vectorgpuField& U,
    const scalar rho
) const
{

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorgpuField& alphaZones = alpha_[zoneI];
        const tensorgpuField& betaZones = beta_[zoneI];

        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();
  
        if(alphaZones.size() == 1)
        {
            tensor alpha = alphaZones.get(0);
            tensor beta = betaZones.get(0);

            thrust::transform
            (
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.begin()
                ),
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.end()
                ),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        U.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.begin()
                ),
                fixeCoeffConstAUFunctor(alpha,beta,rho)
            );
        }
        else
        {
            thrust::transform
            (
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.begin()
                ),
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.end()
                ),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    alphaZones.begin(),
                    betaZones.begin(),
                    thrust::make_permutation_iterator
                    (
                        U.begin(),
                        cells.begin()
                    )
                )),
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.begin()
                ),
                fixeCoeffAUFunctor(rho)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    alphaXYZ_(coeffs_.lookup("alpha")),
    betaXYZ_(coeffs_.lookup("beta")),
    alpha_(cellZoneIDs_.size()),
    beta_(cellZoneIDs_.size())
{
    adjustNegativeResistance(alphaXYZ_);
    adjustNegativeResistance(betaXYZ_);

    calcTranformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::~fixedCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::calcTranformModelData()
{
    if (coordSys_.R().uniform())
    {
        forAll (cellZoneIDs_, zoneI)
        {
            alpha_[zoneI].setSize(1);
            beta_[zoneI].setSize(1);

            tensor a = tensor::zero;
            a.xx() = alphaXYZ_.value().x();
            a.yy() = alphaXYZ_.value().y();
            a.zz() = alphaXYZ_.value().z();
            
            alpha_[zoneI].set(0,coordSys_.R().transformTensor(a));

            tensor b = tensor::zero;
            b.xx() = betaXYZ_.value().x();
            b.yy() = betaXYZ_.value().y();
            b.zz() = betaXYZ_.value().z();
            beta_[zoneI].set(0,coordSys_.R().transformTensor(b));
        }
    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();

            alpha_[zoneI].setSize(cells.size());
            beta_[zoneI].setSize(cells.size());
            
            tensor a = tensor::zero;
            a.xx() = alphaXYZ_.value().x();
            a.yy() = alphaXYZ_.value().y();
            a.zz() = alphaXYZ_.value().z();
            
            alpha_[zoneI].operator=(a);

            tensor b = tensor::zero;
            b.xx() = betaXYZ_.value().x();
            b.yy() = betaXYZ_.value().y();
            b.zz() = betaXYZ_.value().z();
            beta_[zoneI].operator=(b);

            const coordinateRotation& R = coordSys_.R(mesh_, cells);

            alpha_[zoneI] = R.transformTensor(alpha_[zoneI], cells);
            beta_[zoneI] = R.transformTensor(beta_[zoneI], cells);
        }
    }
}


void Foam::porosityModels::fixedCoeff::calcForce
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
    scalar rhoRef = readScalar(coeffs_.lookup("rhoRef"));

    apply(Udiag, Usource, V, U.getField(), rhoRef);

    force = Udiag*U.getField() - Usource;
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorgpuField& U = UEqn.psi().getField();
    const scalargpuField& V = mesh_.V().getField();
    scalargpuField& Udiag = UEqn.diag();
    vectorgpuField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorgpuField& U = UEqn.psi().getField();
    const scalargpuField& V = mesh_.V().getField();
    scalargpuField& Udiag = UEqn.diag();
    vectorgpuField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorgpuField& U = UEqn.psi().getField();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(AU.getField(), U, rho);
}


bool Foam::porosityModels::fixedCoeff::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
