/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "actuationDiskSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

namespace Foam
{
struct actuationDiskSourceFunctor{
	const scalar V;
	const scalar T;
	const tensor E;
	const vector upU;
	actuationDiskSourceFunctor(const scalar _V, const scalar _T, const tensor _E,const vector _upU):V(_V),T(_T),E(_E),upU(_upU){}
	__HOST____DEVICE__
	vector operator () (const vector& Usource, const scalar& Vcell){
		return Usource + (((Vcell/V)*T*E) & upU);
	}
};
}

template<class RhoFieldType>
void Foam::fv::actuationDiskSource::addActuationDiskAxialInertialResistance
(
    vectorgpuField& Usource,
    const labelgpuList& cells,
    const scalargpuField& Vcells,
    const RhoFieldType& rho,
    const vectorgpuField& U
) const
{
    scalar a = 1.0 - Cp_/Ct_;
    vector uniDiskDir = diskDir_/mag(diskDir_);
    tensor E(tensor::zero);
    E.xx() = uniDiskDir.x();
    E.yy() = uniDiskDir.y();
    E.zz() = uniDiskDir.z();

    vector upU = vector(VGREAT, VGREAT, VGREAT);
    scalar upRho = VGREAT;
    if (upstreamCellId_ != -1)
    {
        upU =  U.get(upstreamCellId_);
        upRho = rho.getField().get(upstreamCellId_);
    }
    reduce(upU, minOp<vector>());
    reduce(upRho, minOp<scalar>());

    scalar T = 2.0*upRho*diskArea_*mag(upU)*a*(1 - a);

/*
    forAll(cells, i)
    {
        Usource[cells[i]] += ((Vcells[cells[i]]/V())*T*E) & upU;
    }
*/  
    thrust::transform(thrust::make_permutation_iterator(Usource.begin(),cells.begin()),
                      thrust::make_permutation_iterator(Usource.begin(),cells.end()),
                      thrust::make_permutation_iterator(Vcells.begin(),cells.begin()),
                      thrust::make_permutation_iterator(Usource.begin(),cells.begin()),
                      actuationDiskSourceFunctor(V(),T,E,upU));
}


// ************************************************************************* //
