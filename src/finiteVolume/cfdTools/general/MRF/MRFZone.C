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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "faceSet.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFZone, 0);
    
    struct MRFZoneCorrectVelocityFunctor : public std::unary_function<vector,vector>     
    {    
        const vector omega;
        const vector origin;

        MRFZoneCorrectVelocityFunctor(vector _omega, vector _origin): 
                                     omega(_omega), origin(_origin) {}
        __HOST____DEVICE__
        vector operator () (const vector& x)
        {
            return omega ^ (x - origin);
        }    
    };
	    
    template<bool rhs>    
    struct MRFZoneAddCoriollisFunctor : public std::binary_function<vector,thrust::tuple<scalar,vector>,vector>    
    {
        const vector omega;

        MRFZoneAddCoriollisFunctor(vector _omega): omega(_omega){}

        __HOST____DEVICE__
        vector operator () (const vector& Usource,const thrust::tuple<scalar,vector>& t)
        {
            vector delta = thrust::get<0>(t)*(omega ^ thrust::get<1>(t));
            if(rhs)
                return Usource + delta;
            else
                return Usource - delta;
        }    
    };
    
    template<bool rhs>    
    struct MRFZoneAddRhoCoriollisFunctor : public std::binary_function<vector,thrust::tuple<scalar,scalar,vector>,vector>     
    {
        const vector omega;

        MRFZoneAddRhoCoriollisFunctor(vector _omega): omega(_omega){}

        __HOST____DEVICE__
        vector operator () (const vector& Usource,const thrust::tuple<scalar,scalar,vector>& t)
        {
            vector delta = thrust::get<0>(t)*thrust::get<1>(t)*(omega ^ thrust::get<2>(t));
            if(rhs)
                return Usource + delta;
            else
                return Usource - delta;
        }    
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0:not in zone
    //  1:moving with frame
    //  2:other
    labelList faceType(mesh_.nFaces(), 0);

    // Determine faces in cell zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    if (cellZoneID_ != -1)
    {
        const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
        forAll(cellLabels, i)
        {
            zoneCell[cellLabels[i]] = true;
        }
    }


    label nZoneFaces = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (zoneCell[own[faceI]] || zoneCell[nei[faceI]])
        {
            faceType[faceI] = 1;
            nZoneFaces++;
        }
    }


    labelHashSet excludedPatches(excludedPatchLabels_);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled() || excludedPatches.found(patchI))
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 2;
                    nZoneFaces++;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;

                if (zoneCell[own[faceI]])
                {
                    faceType[faceI] = 1;
                    nZoneFaces++;
                }
            }
        }
    }

    // Now we have for faceType:
    //  0   : face not in cellZone
    //  1   : internal face or normal patch face
    //  2   : coupled patch face or excluded patch face

    // Sort into lists per patch.

    labelList internalFacesTmp(mesh_.nFaces());
    label nInternal = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceType[faceI] == 1)
        {
            internalFacesTmp[nInternal++] = faceI;
        }
    }

    internalFacesTmp.setSize(nInternal);
    internalFaces_ = internalFacesTmp;

    labelList nIncludedFaces(patches.size(), 0);
    labelList nExcludedFaces(patches.size(), 0);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label faceI = pp.start() + patchFacei;

            if (faceType[faceI] == 1)
            {
                nIncludedFaces[patchi]++;
            }
            else if (faceType[faceI] == 2)
            {
                nExcludedFaces[patchi]++;
            }
        }
    }

    includedFaces_.setSize(patches.size());
    excludedFaces_.setSize(patches.size());
    forAll(nIncludedFaces, patchi)
    {
        includedFaces_[patchi].setSize(nIncludedFaces[patchi]);
        excludedFaces_[patchi].setSize(nExcludedFaces[patchi]);
    }
    nIncludedFaces = 0;
    nExcludedFaces = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        
        labelList includedFacesTmp(includedFaces_[patchi].size());
        labelList excludedFacesTmp(excludedFaces_[patchi].size());

        forAll(pp, patchFacei)
        {
            label faceI = pp.start() + patchFacei;

            if (faceType[faceI] == 1)
            {
                includedFacesTmp[nIncludedFaces[patchi]++] = patchFacei;
            }
            else if (faceType[faceI] == 2)
            {
                excludedFacesTmp[nExcludedFaces[patchi]++] = patchFacei;
            }
        }
        
        includedFaces_[patchi] = includedFacesTmp;
        excludedFaces_[patchi] = excludedFacesTmp;
    }

//TODO implement afte added asList to gpuList
/*
    if (debug)
    {
        faceSet internalFaces(mesh_, "internalFaces", internalFaces_);
        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;
        internalFaces.write();

        faceSet MRFFaces(mesh_, "includedFaces", 100);
        forAll(includedFaces_, patchi)
        {
            forAll(includedFaces_[patchi], i)
            {
                label patchFacei = includedFaces_[patchi][i];
                MRFFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << MRFFaces.size()
            << " patch faces in MRF zone to faceSet "
            << MRFFaces.name() << endl;
        MRFFaces.write();

        faceSet excludedFaces(mesh_, "excludedFaces", 100);
        forAll(excludedFaces_, patchi)
        {
            forAll(excludedFaces_[patchi], i)
            {
                label patchFacei = excludedFaces_[patchi][i];
                excludedFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << excludedFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << excludedFaces.name() << endl;
        excludedFaces.write();
    }
*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    mesh_(mesh),
    name_(name),
    coeffs_(dict),
    active_(true),
    cellZoneName_(cellZoneName),
    cellZoneID_(),
    excludedPatchNames_
    (
        coeffs_.lookupOrDefault("nonRotatingPatches", wordList(0))
    ),
    origin_(coeffs_.lookup("origin")),
    axis_(coeffs_.lookup("axis")),
    omega_(DataEntry<scalar>::New("omega", coeffs_))
{
    if (cellZoneName_ == word::null)
    {
        coeffs_.lookup("active") >> active_;
        coeffs_.lookup("cellZone") >> cellZoneName_;
    }

    if (!active_)
    {
        cellZoneID_ = -1;
    }
    else
    {
        cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        axis_ = axis_/mag(axis_);

        excludedPatchLabels_.setSize(excludedPatchNames_.size());

        forAll(excludedPatchNames_, i)
        {
            excludedPatchLabels_[i] =
                patches.findPatchID(excludedPatchNames_[i]);

            if (excludedPatchLabels_[i] == -1)
            {
                FatalErrorIn
                (
                    "MRFZone"
                    "("
                        "const word&, "
                        "const fvMesh&, "
                        "const dictionary&, "
                        "const word&"
                    ")"
                )
                    << "cannot find MRF patch " << excludedPatchNames_[i]
                    << exit(FatalError);
            }
        }

        bool cellZoneFound = (cellZoneID_ != -1);

        if (!cellZoneFound)
        {
            FatalErrorIn
            (
                "MRFZone"
                "("
                    "const word&, "
                    "const fvMesh&, "
                    "const dictionary&, "
                    "const word&"
                ")"
            )
                << "cannot find MRF cellZone " << cellZoneName_
                << exit(FatalError);
        }

        setMRFFaces();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFZone::Omega() const
{
    return omega_->value(mesh_.time().timeOutputValue())*axis_;
}


void Foam::MRFZone::addCoriolis
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelgpuList& cells = mesh_.cellZones()[cellZoneID_].getList();
    vectorgpuField& ddtUc = ddtU.internalField();

    const vector Omega = this->Omega();

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            ddtUc.begin(),
            cells.begin()
        ),
        thrust::make_permutation_iterator
        (
            ddtUc.begin(),
            cells.end()
        ),
        thrust::make_transform_iterator
        (
            thrust::make_permutation_iterator
            (
                ddtUc.begin(),
                cells.begin()
            ),
            crossOperatorSFFunctor<vector,vector,vector>(Omega)
        ),
        thrust::make_permutation_iterator
        (
            ddtUc.begin(),
            cells.begin()
        ),
        addOperatorFunctor<vector,vector,vector>()
    );
}


void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn, const bool rhs) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelgpuList& cells = mesh_.cellZones()[cellZoneID_].getList();
    const scalargpuField& V = mesh_.V().getField();
    vectorgpuField& Usource = UEqn.source();
    const vectorgpuField& U = UEqn.psi().getField();

    const vector Omega = this->Omega();

    if (rhs)
    {
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.begin()
            ),
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator
                (
                    V.begin(),
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
                Usource.begin(),
                cells.begin()
            ),
            MRFZoneAddCoriollisFunctor<true>(Omega)
        );
    }
    else
    {
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.begin()
            ),
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator
                (
                    V.begin(),
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
                Usource.begin(),
                cells.begin()
            ),
            MRFZoneAddCoriollisFunctor<false>(Omega)
        );
    }
}


void Foam::MRFZone::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    const bool rhs
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelgpuList& cells = mesh_.cellZones()[cellZoneID_].getList();
    const scalargpuField& V = mesh_.V().getField();
    vectorgpuField& Usource = UEqn.source();
    const vectorgpuField& U = UEqn.psi().getField();

    const vector Omega = this->Omega();

    if (rhs)
    {
         thrust::transform
         (
             thrust::make_permutation_iterator
             (
                 Usource.begin(),
                 cells.begin()
             ),
             thrust::make_permutation_iterator
             (
                 Usource.begin(),
                 cells.end()
             ),
             thrust::make_zip_iterator(thrust::make_tuple
             (
                 thrust::make_permutation_iterator
                 (
                     V.begin(),
                     cells.begin()
                 ),
                 thrust::make_permutation_iterator
                 (
                     rho.getField().begin(),
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
                 Usource.begin(),
                 cells.begin()
             ),
             MRFZoneAddRhoCoriollisFunctor<true>(Omega)
         );
    }
    else
    {
        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.begin()
            ),
            thrust::make_permutation_iterator
            (
                Usource.begin(),
                cells.end()
            ),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_permutation_iterator
                (
                    V.begin(),
                    cells.begin()
                ),
                thrust::make_permutation_iterator
                (
                    rho.getField().begin(),
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
                Usource.begin(),
                cells.begin()
            ),
            MRFZoneAddRhoCoriollisFunctor<false>(Omega)
        );
    }
}


void Foam::MRFZone::makeRelative(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();

    const labelgpuList& cells = mesh_.cellZones()[cellZoneID_].getList();

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        thrust::make_transform_iterator
        (
            thrust::make_permutation_iterator
            (
                C.getField().begin(),
                cells.begin()
            ),
            MRFZoneCorrectVelocityFunctor(Omega,origin_)
        ),
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        subtractOperatorFunctor<vector,vector,vector>()
    );

    // Included patches
    forAll(includedFaces_, patchi)
    {
        const labelgpuList& faces = includedFaces_[patchi];

        thrust::fill
        (
            thrust::make_permutation_iterator
            (
                U.boundaryField()[patchi].begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                U.boundaryField()[patchi].begin(),
                faces.end()
            ),
            vector::zero
        );
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        const vectorgpuField& patchC = mesh_.C().boundaryField()[patchi];

        const vectorgpuField& pfld = U.boundaryField()[patchi];
        const labelgpuList& faces = excludedFaces_[patchi];

        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.end()
            ),
            thrust::make_transform_iterator
            (
                thrust::make_permutation_iterator
                (
                    patchC.begin(),
                    faces.begin()
                ),
                MRFZoneCorrectVelocityFunctor(Omega,origin_)
            ),
            thrust::make_permutation_iterator
            (
                U.boundaryField()[patchi].begin(),
                faces.begin()
            ),
            subtractOperatorFunctor<vector,vector,vector>()
        );
    }
}


void Foam::MRFZone::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::MRFZone::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::MRFZone::makeAbsolute(volVectorField& U) const
{
    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();

    const labelgpuList& cells = mesh_.cellZones()[cellZoneID_].getList();

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        thrust::make_transform_iterator
        (
            thrust::make_permutation_iterator
            (
               C.getField().begin(),
               cells.begin()
            ),
            MRFZoneCorrectVelocityFunctor(Omega,origin_)
        ),
        thrust::make_permutation_iterator
        (
            U.getField().begin(),
            cells.begin()
        ),
        addOperatorFunctor<vector,vector,vector>()
    );

    // Included patches
    forAll(includedFaces_, patchi)
    {
        const vectorgpuField& patchC = C.boundaryField()[patchi];

        vectorgpuField& pfld = U.boundaryField()[patchi];
        const labelgpuList& faces = includedFaces_[patchi];

        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                patchC.begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                patchC.begin(),
                faces.end()
            ),
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.begin()
            ),
            MRFZoneCorrectVelocityFunctor(Omega,origin_)
        );
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        const vectorgpuField& patchC = mesh_.C().boundaryField()[patchi];

        const vectorgpuField& pfld = U.boundaryField()[patchi];
        const labelgpuList& faces = excludedFaces_[patchi];

        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.end()
            ),
            thrust::make_transform_iterator
            (
                thrust::make_permutation_iterator
                (
                    patchC.begin(),
                    faces.begin()
                ),
                MRFZoneCorrectVelocityFunctor(Omega,origin_)
            ),
            thrust::make_permutation_iterator
            (
                U.boundaryField()[patchi].begin(),
                faces.begin()
            ),
            addOperatorFunctor<vector,vector,vector>()
        );
              
    }
}


void Foam::MRFZone::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    const vector Omega = this->Omega();


    // Included patches
    forAll(includedFaces_, patchi)
    {
        const vectorgpuField& patchC = mesh_.Cf().boundaryField()[patchi];

        vectorgpuField pfld(U.boundaryField()[patchi]);
        const labelgpuList& faces = includedFaces_[patchi];

        thrust::transform
        (
            thrust::make_permutation_iterator
            (
                patchC.begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                patchC.begin(),
                faces.end()
            ),
            thrust::make_permutation_iterator
            (
                pfld.begin(),
                faces.begin()
            ),
            MRFZoneCorrectVelocityFunctor(Omega,origin_)
        );

        U.boundaryField()[patchi] == pfld;
        
    }
}


void Foam::MRFZone::writeData(Ostream& os) const
{
    os  << nl;
    os.write(name_) << nl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("active") << active_ << token::END_STATEMENT << nl;
    os.writeKeyword("cellZone") << cellZoneName_ << token::END_STATEMENT << nl;
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    omega_->writeData(os);

    if (excludedPatchNames_.size())
    {
        os.writeKeyword("nonRotatingPatches") << excludedPatchNames_
            << token::END_STATEMENT << nl;
    }

    os  << decrIndent << token::END_BLOCK << nl;
}


bool Foam::MRFZone::read(const dictionary& dict)
{
    coeffs_ = dict;

    active_ = readBool(coeffs_.lookup("active"));
    coeffs_.lookup("cellZone") >> cellZoneName_;
    cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

    return true;
}


// ************************************************************************* //
