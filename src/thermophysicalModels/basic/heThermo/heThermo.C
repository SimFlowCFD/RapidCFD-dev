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

#include "heThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace Foam
{
	template<class Mixture>
	struct heThermoHEFunctor{
		const Mixture mixture;
		heThermoHEFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.HE(p,T);
		}
	};
	
	template<class Mixture>
	struct heThermoCpFunctor{
		const Mixture mixture;
		heThermoCpFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.Cp(p,T);
		}
	};	
	
	template<class Mixture>
	struct heThermoCvFunctor{
		const Mixture mixture;
		heThermoCvFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.Cv(p,T);
		}
	};
	
    template<class Mixture>
	struct heThermoGammaFunctor{
		const Mixture mixture;
		heThermoGammaFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.gamma(p,T);
		}
	};
	
	template<class Mixture>
	struct heThermoCpvFunctor{
		const Mixture mixture;
		heThermoCpvFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.Cpv(p,T);
		}
	};
	
	template<class Mixture>
	struct heThermoCpByCpvFunctor{
		const Mixture mixture;
		heThermoCpByCpvFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const scalar& p, const scalar& T){
			return mixture.cpBycpv(p,T);
		}
	};
	
	template<class Mixture>
	struct heThermoTHEFunctor{
		const Mixture mixture;
		heThermoTHEFunctor(const Mixture _mixture): mixture(_mixture) {}
		__HOST____DEVICE__
		scalar operator () (const thrust::tuple<scalar,scalar,scalar>& t){
			const scalar h = thrust::get<0>(t);
			const scalar p = thrust::get<1>(t);
			const scalar T = thrust::get<2>(t);
			return mixture.THE(h,p,T);
		}
	};
}

template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::
heBoundaryCorrection(volScalarField& h)
{
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::init()
{
    scalargpuField& heCells = he_.internalField();
    const scalargpuField& pCells = this->p_.internalField();
    const scalargpuField& TCells = this->T_.internalField();
/*
    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TCells[celli]);
    }
*/
    thrust::transform(pCells.begin(),pCells.end(),TCells.begin(),heCells.begin(),
                      heThermoHEFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(he_.boundaryField(), patchi)
    {
        he_.boundaryField()[patchi] == he
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        );
    }

    this->heBoundaryCorrection(he_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BasicThermo(mesh, phaseName),
    MixtureType(*this, mesh),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init();
}


template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::~heThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> the
    (
        new volScalarField
        (
            IOobject
            (
                "he",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& he = the();
    scalargpuField& heCells = he.internalField();
    const scalargpuField& pCells = p.internalField();
    const scalargpuField& TCells = T.internalField();
/*
    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TCells[celli]);
    }
*/
    thrust::transform(pCells.begin(),pCells.end(),TCells.begin(),heCells.begin(),
                      heThermoHEFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(he.boundaryField(), patchi)
    {
        scalargpuField& hep = he.boundaryField()[patchi];
        const scalargpuField& pp = p.boundaryField()[patchi];
        const scalargpuField& Tp = T.boundaryField()[patchi];
/*
        forAll(hep, facei)
        {
            hep[facei] =
                this->patchFaceMixture(patchi, facei).HE(pp[facei], Tp[facei]);
        }
*/
        thrust::transform(pp.begin(),pp.end(),Tp.begin(),hep.begin(),
                      heThermoHEFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi, 0)));
    }

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalargpuField& p,
    const scalargpuField& T,
    const labelgpuList& cells
) const
{
    tmp<scalargpuField> the(new scalargpuField(T.size()));
    scalargpuField& he = the();
/*
    forAll(T, celli)
    {
        he[celli] = this->cellMixture(cells[celli]).HE(p[celli], T[celli]);
    }
*/
    thrust::transform(p.begin(),p.end(),T.begin(),he.begin(),
                      heThermoHEFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> the(new scalargpuField(T.size()));
    scalargpuField& he = the();
/*
    forAll(T, facei)
    {
        he[facei] =
            this->patchFaceMixture(patchi, facei).HE(p[facei], T[facei]);
    }
*/
    thrust::transform(p.begin(),p.end(),T.begin(),he.begin(),
                      heThermoHEFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi, 0)));

    return the;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            he_.dimensions()
        )
    );

    volScalarField& hcf = thc();
    scalargpuField& hcCells = hcf.internalField();
    
    scalar hc = this->cellMixture(0).Hc();
    hcCells = hc;
/*
    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }
*/
    forAll(hcf.boundaryField(), patchi)
    {
        scalargpuField& hcp = hcf.boundaryField()[patchi];
        hcp = hc;
/*
        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
*/
    }

    return thc;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::Cp
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> tCp(new scalargpuField(T.size()));
    scalargpuField& cp = tCp();
/*
    forAll(T, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp(p[facei], T[facei]);
    }
*/
    thrust::transform(p.begin(),p.end(),T.begin(),cp.begin(),
                      heThermoCpFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp();
/*
    forAll(this->T_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp(this->p_[celli], this->T_[celli]);
    }
*/    
    thrust::transform(this->p_.getField().begin(),this->p_.getField().end(),this->T_.getField().begin(),cp.getField().begin(),
                      heThermoCpFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));
    

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];
/*
        forAll(pT, facei)
        {
            pCp[facei] =
                this->patchFaceMixture(patchi, facei).Cp(pp[facei], pT[facei]);
        }
*/
        thrust::transform(pp.begin(),pp.end(),pT.begin(),pCp.begin(),
                      heThermoCpFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField>
Foam::heThermo<BasicThermo, MixtureType>::Cv
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> tCv(new scalargpuField(T.size()));
    scalargpuField& cv = tCv();
/*
    forAll(T, facei)
    {
        cv[facei] =
            this->patchFaceMixture(patchi, facei).Cv(p[facei], T[facei]);
    }
*/
    thrust::transform(p.begin(),p.end(),T.begin(),cv.begin(),
                      heThermoCpFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv();
/*
    forAll(this->T_, celli)
    {
        cv[celli] =
            this->cellMixture(celli).Cv(this->p_[celli], this->T_[celli]);
    }
*/
    thrust::transform(this->p_.getField().begin(),this->p_.getField().end(),this->T_.getField().begin(),cv.getField().begin(),
                      heThermoCvFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(this->T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        );
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::gamma
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> tgamma(new scalargpuField(T.size()));
    scalargpuField& cpv = tgamma();
/*
    forAll(T, facei)
    {
        cpv[facei] =
            this->patchFaceMixture(patchi, facei).gamma(p[facei], T[facei]);
    }
*/     
    thrust::transform(p.begin(),p.end(),T.begin(),cpv.begin(),
                      heThermoGammaFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tgamma
    (
        new volScalarField
        (
            IOobject
            (
                "gamma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& cpv = tgamma();
/*
    forAll(this->T_, celli)
    {
        cpv[celli] =
            this->cellMixture(celli).gamma(this->p_[celli], this->T_[celli]);
    }
*/
    thrust::transform(this->p_.getField().begin(),this->p_.getField().end(),this->T_.getField().begin(),cpv.getField().begin(),
                      heThermoGammaFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = cpv.boundaryField()[patchi];
/*
        forAll(pT, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                pp[facei],
                pT[facei]
            );
        }
*/
        thrust::transform(pp.begin(),pp.end(),pT.begin(),pgamma.begin(),
                      heThermoGammaFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));
    }

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::Cpv
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> tCpv(new scalargpuField(T.size()));
    scalargpuField& cpv = tCpv();
/*
    forAll(T, facei)
    {
        cpv[facei] =
            this->patchFaceMixture(patchi, facei).Cpv(p[facei], T[facei]);
    }
*/
    thrust::transform(p.begin(),p.end(),T.begin(),cpv.begin(),
                      heThermoCpvFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));

    return tCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cpv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCpv
    (
        new volScalarField
        (
            IOobject
            (
                "Cpv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cpv = tCpv();
/*
    forAll(this->T_, celli)
    {
        cpv[celli] =
            this->cellMixture(celli).Cpv(this->p_[celli], this->T_[celli]);
    }
*/
    thrust::transform(this->p_.getField().begin(),this->p_.getField().end(),this->T_.getField().begin(),cpv.getField().begin(),
                      heThermoCpvFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpv = cpv.boundaryField()[patchi];
/*
        forAll(pT, facei)
        {
            pCpv[facei] =
                this->patchFaceMixture(patchi, facei).Cpv(pp[facei], pT[facei]);
        }
*/
        thrust::transform(pp.begin(),pp.end(),pT.begin(),pCpv.begin(),
                      heThermoCpvFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));
        
    }

    return tCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::CpByCpv
(
    const scalargpuField& p,
    const scalargpuField& T,
    const label patchi
) const
{
    tmp<scalargpuField> tCpByCpv(new scalargpuField(T.size()));
    scalargpuField& cpByCpv = tCpByCpv();
/*
    forAll(T, facei)
    {
        cpByCpv[facei] =
            this->patchFaceMixture(patchi, facei).cpBycpv(p[facei], T[facei]);
    }
*/

    thrust::transform(p.begin(),p.end(),T.begin(),cpByCpv.begin(),
                      heThermoCpByCpvFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));    
    

    return tCpByCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::CpByCpv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCpByCpv
    (
        new volScalarField
        (
            IOobject
            (
                "CpByCpv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& cpByCpv = tCpByCpv();
/*
    forAll(this->T_, celli)
    {
        cpByCpv[celli] = this->cellMixture(celli).cpBycpv
        (
            this->p_[celli],
            this->T_[celli]
        );
    }
*/
    thrust::transform(this->p_.getField().begin(),this->p_.getField().end(),this->T_.getField().begin(),cpByCpv.getField().begin(),
                      heThermoCpByCpvFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCpByCpv = cpByCpv.boundaryField()[patchi];
/*
        forAll(pT, facei)
        {
            pCpByCpv[facei] = this->patchFaceMixture(patchi, facei).cpBycpv
            (
                pp[facei],
                pT[facei]
            );
        }
*/
        thrust::transform(pp.begin(),pp.end(),pT.begin(),pCpByCpv.begin(),
                      heThermoCpvFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));
    }

    return tCpByCpv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalargpuField& h,
    const scalargpuField& p,
    const scalargpuField& T0,
    const labelgpuList& cells
) const
{
    tmp<scalargpuField> tT(new scalargpuField(h.size()));
    scalargpuField& T = tT();
/*
    forAll(h, celli)
    {
        T[celli] =
            this->cellMixture(cells[celli]).THE(h[celli], p[celli], T0[celli]);
    }
*/
    thrust::transform(thrust::make_zip_iterator(
                              thrust::make_tuple(
                                      h.begin(),
                                      p.begin(),
                                      T0.begin()
                              )),
                      thrust::make_zip_iterator(
                              thrust::make_tuple(
                                      h.end(),
                                      p.end(),
                                      T0.end()
                              )),     
                      T.begin(),
                      heThermoTHEFunctor< typename MixtureType::thermoType >(this->cellMixture(0)));

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalargpuField& h,
    const scalargpuField& p,
    const scalargpuField& T0,
    const label patchi
) const
{

    tmp<scalargpuField> tT(new scalargpuField(h.size()));
    scalargpuField& T = tT();
    /*
    forAll(h, facei)
    {
        T[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).THE(h[facei], p[facei], T0[facei]);
    }
    */
    
    thrust::transform(thrust::make_zip_iterator(
                              thrust::make_tuple(
                                      h.begin(),
                                      p.begin(),
                                      T0.begin()
                              )),
                      thrust::make_zip_iterator(
                              thrust::make_tuple(
                                      h.end(),
                                      p.end(),
                                      T0.end()
                              )),     
                      T.begin(),
                      heThermoTHEFunctor< typename MixtureType::thermoType >(this->patchFaceMixture(patchi,0)));

    return tT;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappa() const
{
    tmp<Foam::volScalarField> kappa(Cp()*this->alpha_);
    kappa().rename("kappa");
    return kappa;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField> Foam::heThermo<BasicThermo, MixtureType>::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(Cp()*alphaEff(alphat));
    kappaEff().rename("kappaEff");
    return kappaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalargpuField& alphat,
    const label patchi
) const
{
    return
        Cp
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        )*alphaEff(alphat, patchi);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpv()*(this->alpha_ + alphat));
    alphaEff().rename("alphaEff");
    return alphaEff;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalargpuField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalargpuField& alphat,
    const label patchi
) const
{
    return
    this->CpByCpv
    (
        this->p_.boundaryField()[patchi],
        this->T_.boundaryField()[patchi],
        patchi
    )
   *(
        this->alpha_.boundaryField()[patchi]
      + alphat
    );
}


template<class BasicThermo, class MixtureType>
bool Foam::heThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
