/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    struct powerLawUdiagFunctor
    {
        const scalar C0;
        const scalar C1m1b2;
        powerLawUdiagFunctor(scalar _C0,scalar _C1m1b2): 
                            C0(_C0), C1m1b2(_C1m1b2){}
        __HOST____DEVICE__                                          //V       rho     U
        scalar operator () (const scalar& Udiag,const thrust::tuple<scalar, scalar, vector>& t)
        {
            return Udiag + thrust::get<0>(t)*thrust::get<1>(t)*C0*pow(magSqr(thrust::get<2>(t)),C1m1b2);
        }
    };
	
    struct powerLawAUFunctor
    {
        const scalar C0;
        const scalar C1m1b2;
        const tensor I;
        powerLawAUFunctor(scalar _C0,scalar _C1m1b2): 
                      C0(_C0), C1m1b2(_C1m1b2),I(tensor::zero){}
        __HOST____DEVICE__                                        //rho     U
        tensor operator () (const tensor& AU,const thrust::tuple<scalar, vector>& t)
        {
            return AU + I*(thrust::get<0>(t)*C0*pow(magSqr(thrust::get<1>(t)),C1m1b2));
        }
    };
	
}

template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    scalargpuField& Udiag,
    const scalargpuField& V,
    const RhoFieldType& rho,
    const vectorgpuField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();
  
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                Udiag.begin(),
                cells.begin()
            ),
            thrust::make_permutation_iterator
            (
                Udiag.begin(),
                cells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator
                (
                V.begin(),cells.begin()
                ),
                thrust::make_permutation_iterator
                (
                rho.begin(),cells.begin()
                ),
                thrust::make_permutation_iterator
                (
                U.begin(),cells.begin()
                )
            )),
            thrust::make_permutation_iterator
            (
                Udiag.begin(),
                cells.begin()
            ),
            powerLawUdiagFunctor(C0,C1m1b2)
        );
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    tensorgpuField& AU,
    const RhoFieldType& rho,
    const vectorgpuField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cellZoneIDs_, zoneI)
    {
        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();
     
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
                    rho.begin(),
                    cells.begin()
                ),
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
            powerLawAUFunctor(C0,C1m1b2)
        );
    }
}


// ************************************************************************* //
