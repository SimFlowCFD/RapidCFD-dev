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

#include "hePsiThermo.H"

namespace Foam
{
	template<class MixtureType>
	struct hePsiThermoCalculateFunctor{
		const MixtureType mixture;
		hePsiThermoCalculateFunctor(const MixtureType _mixture): mixture(_mixture){}
		__HOST____DEVICE__
		thrust::tuple<scalar,scalar,scalar,scalar>
		operator ()(const scalar& h, const thrust::tuple<scalar,scalar>& t){
			scalar p = thrust::get<0>(t);
			scalar T = mixture.THE(h,p,thrust::get<1>(t));
			
			return thrust::make_tuple(T,
			                          mixture.psi(p,T),
			                          mixture.mu(p,T),
			                          mixture.alphah(p,T)
			                         );
		}
	};
	
	template<class MixtureType>
	struct hePsiThermoHECalculateFunctor{
		const MixtureType mixture;
		hePsiThermoHECalculateFunctor(const MixtureType _mixture): mixture(_mixture){}
		__HOST____DEVICE__
		thrust::tuple<scalar,scalar,scalar,scalar>
		operator ()(const scalar& p, const scalar& T){
			
			return thrust::make_tuple(mixture.HE(p,T),
			                          mixture.psi(p,T),
			                          mixture.mu(p,T),
			                          mixture.alphah(p,T)
			                         );
		}
	};
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalargpuField& hCells = this->he_.internalField();
    const scalargpuField& pCells = this->p_.internalField();

    scalargpuField& TCells = this->T_.internalField();
    scalargpuField& psiCells = this->psi_.internalField();
    scalargpuField& muCells = this->mu_.internalField();
    scalargpuField& alphaCells = this->alpha_.internalField();   

/*
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }
*/
    thrust::transform(hCells.begin(),hCells.end(),
                      thrust::make_zip_iterator(thrust::make_tuple( pCells.begin(),
                                                                    TCells.begin())),
                      thrust::make_zip_iterator(thrust::make_tuple(TCells.begin(),
                                                                   psiCells.begin(),
                                                                   muCells.begin(),
                                                                   alphaCells.begin()
                                                                   )),
                      hePsiThermoCalculateFunctor<typename MixtureType::thermoType>(this->cellMixture(0)));
                    

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = this->he_.boundaryField()[patchi];

        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
			/*
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
            */
            thrust::transform(pp.begin(),pp.end(),
                              pT.begin(),
					  thrust::make_zip_iterator(thrust::make_tuple(ph.begin(),
																   ppsi.begin(),
																   pmu.begin(),
																   palpha.begin()
																   )),
					  hePsiThermoHECalculateFunctor<typename MixtureType::thermoType>(this->patchFaceMixture(patchi, 0)));

        }
        else
        {
			/*
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.THE(ph[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
            */
            
			thrust::transform(ph.begin(),ph.end(),
					  thrust::make_zip_iterator(thrust::make_tuple( pp.begin(),
																	pT.begin())),
					  thrust::make_zip_iterator(thrust::make_tuple(ph.begin(),
																   ppsi.begin(),
																   pmu.begin(),
																   palpha.begin()
																   )),
					  hePsiThermoCalculateFunctor<typename MixtureType::thermoType>(this->patchFaceMixture(patchi, 0)));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::hePsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::~hePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering hePsiThermo<BasicPsiThermo, MixtureType>::correct()"
            << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting hePsiThermo<BasicPsiThermo, MixtureType>::correct()"
            << endl;
    }
}


// ************************************************************************* //
