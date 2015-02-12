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

#include "heSolidThermo.H"
#include "volFields.H"

namespace Foam
{
	template<class MixtureType>
	struct heSolidThermoCalculateFunctor{
		const MixtureType mixture;
		const MixtureType volMixture;
		heSolidThermoCalculateFunctor(const MixtureType _mixture,const MixtureType _volMixture): mixture(_mixture),volMixture(_volMixture){}
		__HOST____DEVICE__
		thrust::tuple<scalar,scalar,scalar>
		operator ()(const scalar& h, const thrust::tuple<scalar,scalar>& t){
			scalar p = thrust::get<0>(t);
			scalar T = mixture.THE(h,p,thrust::get<1>(t));
			
			return thrust::make_tuple(T,
			                          volMixture.rho(p,T),
			                          volMixture.kappa(p,T)/mixture.Cpv(p,T)
			                         );
		}
	};
	
	template<class MixtureType>
	struct heSolidThermoHECalculateFunctor{
		const MixtureType mixture;
		const MixtureType volMixture;
		heSolidThermoHECalculateFunctor(const MixtureType _mixture,const MixtureType _volMixture): mixture(_mixture),volMixture(_volMixture){}
		__HOST____DEVICE__
		thrust::tuple<scalar,scalar,scalar>
		operator ()(const scalar& p, const scalar& T){
			
			return thrust::make_tuple(mixture.HE(p,T),
			                          volMixture.rho(p,T),
			                          volMixture.kappa(p,T)/mixture.Cpv(p,T)
			                         );
		}
	};
	
	template<class Mixture>
	struct heSolidThermoKappaFunctor{
		const Mixture mixture;
		heSolidThermoKappaFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		vector operator () (const scalar& p, const scalar& T){
			return mixture.Kappa(p,T);
		}
	};
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    scalargpuField& TCells = this->T_.internalField();

    const scalargpuField& hCells = this->he_.internalField();
    const scalargpuField& pCells = this->p_.internalField();
    scalargpuField& rhoCells = this->rho_.internalField();
    scalargpuField& alphaCells = this->alpha_.internalField();
/*
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        const typename MixtureType::thermoType& volMixture_ =
            this->cellVolMixture(pCells[celli], TCells[celli], celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        rhoCells[celli] = volMixture_.rho(pCells[celli], TCells[celli]);

        alphaCells[celli] =
            volMixture_.kappa(pCells[celli], TCells[celli])
            /
            mixture_.Cpv(pCells[celli], TCells[celli]);
    }
*/    
        thrust::transform(hCells.begin(),hCells.end(),
                      thrust::make_zip_iterator(thrust::make_tuple( pCells.begin(),
                                                                    TCells.begin())),
                      thrust::make_zip_iterator(thrust::make_tuple(TCells.begin(),
                                                                   rhoCells.begin(),
                                                                   alphaCells.begin()
                                                                   )),
                      heSolidThermoCalculateFunctor<typename MixtureType::thermoType>(this->cellMixture(0),
                                                                                      this->cellVolMixture(0, 0, 0)));

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        fvPatchScalarField& ph = this->he_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            /*
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );


                ph[facei] = mixture_.HE(pp[facei], pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    volMixture_.kappa(pp[facei], pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
            }
            */
            thrust::transform(pp.begin(),pp.end(),
                              pT.begin(),
					  thrust::make_zip_iterator(thrust::make_tuple(ph.begin(),
																   prho.begin(),
																   palpha.begin()
																   )),
					  heSolidThermoHECalculateFunctor<typename MixtureType::thermoType>(this->patchFaceMixture(patchi, 0),
					                                                                  this->patchFaceVolMixture(0,0,patchi, 0)));
		    
        }
        else
        {
			/*
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );

                pT[facei] = mixture_.THE(ph[facei], pp[facei] ,pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    volMixture_.kappa(pp[facei], pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
            }
            */
            thrust::transform(ph.begin(),ph.end(),
					  thrust::make_zip_iterator(thrust::make_tuple( pp.begin(),
																	pT.begin())),
					  thrust::make_zip_iterator(thrust::make_tuple(pT.begin(),
																   prho.begin(),
																   palpha.begin()
																   )),
					  heSolidThermoCalculateFunctor<typename MixtureType::thermoType>(this->patchFaceMixture(patchi, 0),
					                                                                  this->patchFaceVolMixture(0,0,patchi, 0)));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh, dict, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::~heSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heSolidThermo<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting heSolidThermo<MixtureType>::correct()" << endl;
    }
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volVectorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volVectorField> tKappa
    (
        new volVectorField
        (
            IOobject
            (
                "Kappa",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimTime/dimLength/dimTemperature
        )
    );

    volVectorField& Kappa = tKappa();
    vectorgpuField& KappaCells = Kappa.internalField();
    const scalargpuField& TCells = this->T_.internalField();
    const scalargpuField& pCells = this->p_.internalField();
/*
    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).Kappa(pCells[celli], TCells[celli]);
    }
    * */
    thrust::transform(pCells.begin(),pCells.end(),TCells.begin(),KappaCells.begin(),
                      heSolidThermoKappaFunctor<typename MixtureType::thermoType>(this->cellVolMixture(0, 0, 0)));

    forAll(Kappa.boundaryField(), patchi)
    {
        vectorgpuField& Kappap = Kappa.boundaryField()[patchi];
        const scalargpuField& pT = this->T_.boundaryField()[patchi];
        const scalargpuField& pp = this->p_.boundaryField()[patchi];
/*
        forAll(Kappap, facei)
        {
            Kappap[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
                    patchi,
                    facei
                ).Kappa(pp[facei], pT[facei]);
        }
 */       
        thrust::transform(pp.begin(),pp.end(),pT.begin(),Kappap.begin(),
                      heSolidThermoKappaFunctor<typename MixtureType::thermoType>(this->patchFaceVolMixture(0, 0,patchi, 0)));
    }

    return tKappa;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::vectorgpuField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa
(
    const label patchi
) const
{
    const scalargpuField& pp = this->p_.boundaryField()[patchi];
    const scalargpuField& Tp = this->T_.boundaryField()[patchi];
    tmp<vectorgpuField> tKappa(new vectorgpuField(pp.size()));

    vectorgpuField& Kappap = tKappa();
/*
    forAll(Tp, facei)
    {
        Kappap[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).Kappa(pp[facei], Tp[facei]);
    }
*/
    thrust::transform(pp.begin(),pp.end(),Tp.begin(),Kappap.begin(),
                      heSolidThermoKappaFunctor<typename MixtureType::thermoType>(this->patchFaceVolMixture(0, 0,patchi, 0)));

    return tKappa;
}


// ************************************************************************* //
