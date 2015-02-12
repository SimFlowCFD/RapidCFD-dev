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
__HOST____DEVICE__
inline tensor calcCd(const tensor& d, const tensor& f, const scalar& mu, const scalar& rho, const vector& U)
{
    return mu*d + (rho*mag(U))*f;
}

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
inline thrust::tuple<scalar,vector> DFapply(const tensor& d, const tensor& f,const scalar& mu,const scalar& rho,const vector& U,const scalar& V)
{
    const tensor Cd = calcCd(d,f,mu,rho,U);
    const scalar isoCd = tr(Cd);
	
    scalar Udiag = calcUdiag(V,isoCd);
    vector Usource = calcUsource(V,isoCd,Cd,U);
    
    return thrust::make_tuple(Udiag,Usource);
}

struct DarcyForchheimerFunctor
{
    __HOST____DEVICE__
    thrust::tuple<scalar,vector> operator ()(const thrust::tuple<scalar,vector>& in,const thrust::tuple<tensor,tensor,scalar,scalar,vector,scalar>& t)
    {
        thrust::tuple<scalar,vector>  out =  
            DFapply
            (
                thrust::get<0>(t),
                thrust::get<1>(t),
                thrust::get<2>(t),
                thrust::get<3>(t),
                thrust::get<4>(t),
                thrust::get<5>(t)
            );

         return thrust::make_tuple
             (
                 thrust::get<0>(in)+thrust::get<0>(out),
                 thrust::get<1>(in)-thrust::get<1>(out)
             );
    }
};

struct DarcyForchheimerConstFunctor
{
    const tensor d;
    const tensor f;

    DarcyForchheimerConstFunctor(tensor d_, tensor f_): d(d_), f(f_) {}

    __HOST____DEVICE__
    thrust::tuple<scalar,vector> operator ()(const thrust::tuple<scalar,vector>& in,const thrust::tuple<scalar,scalar,vector,scalar>& t)
    {
        thrust::tuple<scalar,vector>  out =  
            DFapply
            (
                d,
                f,
                thrust::get<0>(t),
                thrust::get<1>(t),
                thrust::get<2>(t),
                thrust::get<3>(t)
            );

        return thrust::make_tuple
            (
                thrust::get<0>(in)+thrust::get<0>(out),
                thrust::get<1>(in)-thrust::get<1>(out)
            );
    }
};

struct DarcyForchheimerAUFunctor
{
    __HOST____DEVICE__
    tensor operator ()(const tensor& AU,const thrust::tuple<tensor,tensor,scalar,scalar,vector>& t)
    {
        return AU + calcCd
            (
                thrust::get<0>(t),
                thrust::get<1>(t),
                thrust::get<2>(t),
                thrust::get<3>(t),
                thrust::get<4>(t)
           );
    }
};

struct DarcyForchheimerConstAUFunctor
{
    const tensor d;
    const tensor f;

    DarcyForchheimerConstAUFunctor(tensor d_, tensor f_): d(d_), f(f_) {}
    __HOST____DEVICE__
    tensor operator ()(const tensor& AU,const thrust::tuple<scalar,scalar,vector>& t)
    {
        return AU + calcCd
            (
                d,
                f,
                thrust::get<0>(t),
                thrust::get<1>(t),
                thrust::get<2>(t)
            );
    }
};
}
template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    scalargpuField& Udiag,
    vectorgpuField& Usource,
    const scalargpuField& V,
    const RhoFieldType& rho,
    const scalargpuField& mu,
    const vectorgpuField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorgpuField& dZones = D_[zoneI];
        const tensorgpuField& fZones = F_[zoneI];

        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();

        if(dZones.size() == 1)
        {
            tensor d = dZones.get(0);
            tensor f = fZones.get(0);

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
                        mu.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        rho.begin(),
                        cells.begin()
                    ),
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
                DarcyForchheimerConstFunctor(d,f)
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
                    dZones.begin(),
                    fZones.begin(),
                    thrust::make_permutation_iterator
                    (
                        mu.begin(),
                        cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        rho.begin(),
                        cells.begin()
                    ),
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
                DarcyForchheimerFunctor()
            );
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    tensorgpuField& AU,
    const RhoFieldType& rho,
    const scalargpuField& mu,
    const vectorgpuField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorgpuField& dZones = D_[zoneI];
        const tensorgpuField& fZones = F_[zoneI];

        const labelgpuList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]].getList();

        
        if(dZones.size() == 1)
        {
            tensor d = dZones.get(0);
            tensor f = fZones.get(0);

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
                       mu.begin(),
                       cells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                       rho.begin(),
                       cells.begin()
                    ),
                    thrust::make_permutation_iterator(
                       U.begin(),
                       cells.begin())
                )),
                thrust::make_permutation_iterator
                (
                    AU.begin(),
                    cells.begin()
                ),
                DarcyForchheimerConstAUFunctor(d,f)
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
                   dZones.begin(),
                   fZones.begin(),
                   thrust::make_permutation_iterator
                   (
                       mu.begin(),
                       cells.begin()
                   ),
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
               DarcyForchheimerAUFunctor()
           );
        }
        
    }
}


// ************************************************************************* //
