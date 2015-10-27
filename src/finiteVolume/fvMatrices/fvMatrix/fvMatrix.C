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

#include "volFields.H"
#include "surfaceFields.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "UIndirectList.H"
#include "fvMatrixCache.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    template<class Type,bool add>
    struct fvMatrixPatchAddFunctor
    {
        const Type* issf;
        const label* neiStart;
        const label* losort;

        fvMatrixPatchAddFunctor
        (
            const Type* _issf,
            const label* _neiStart,
            const label* _losort
        ):
             issf(_issf),
             neiStart(_neiStart),
             losort(_losort)
        {}

        __host__ __device__
        Type operator()(const Type& d, const label& id)
        {
            Type out = d;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                if(add)
                    out += issf[face];
                else
                    out -= issf[face];
            }

            return out;
        }
    };
}

template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelgpuList& addr,
    const labelgpuList& sort,
    const labelgpuList& sortStart,
    const gpuField<Type2>& pf,
    gpuField<Type2>& intf
) const
{
    if (sort.size() != pf.size())
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const labelUList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.begin()
        ),
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.end()
        ),
        thrust::make_counting_iterator(0),
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.begin()
        ),
        fvMatrixPatchAddFunctor<Type2,true>
        (
            pf.data(),
            sortStart.data(),
            sort.data()
        )
    );
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::addToInternalField
(
    const labelgpuList& addr,
    const labelgpuList& sort,
    const labelgpuList& sortStart,
    const tmp<gpuField<Type2> >& tpf,
    gpuField<Type2>& intf
) const
{
    addToInternalField(addr,sort,sortStart,tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelgpuList& addr,
    const labelgpuList& sort,
    const labelgpuList& sortStart,
    const gpuField<Type2>& pf,
    gpuField<Type2>& intf
) const
{
    if (sort.size() != pf.size())
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const labelUList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.begin()
        ),
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.end()
        ),
        thrust::make_counting_iterator(0),
        thrust::make_permutation_iterator
        (
            intf.begin(),
            addr.begin()
        ),
        fvMatrixPatchAddFunctor<Type,false>
        (
            pf.data(),
            sortStart.data(),
            sort.data()
        )
    );
}


template<class Type>
template<class Type2>
void Foam::fvMatrix<Type>::subtractFromInternalField
(
    const labelgpuList& addr,
    const labelgpuList& sort,
    const labelgpuList& sortStart,
    const tmp<gpuField<Type2> >& tpf,
    gpuField<Type2>& intf
) const
{
    subtractFromInternalField(addr,sort,sortStart,tpf(), intf);
    tpf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::addBoundaryDiag
(
    scalargpuField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchSortCells(patchI),
            lduAddr().patchSortAddr(patchI),
            lduAddr().patchSortStartAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::addCmptAvBoundaryDiag(scalargpuField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchSortCells(patchI),
            lduAddr().patchSortAddr(patchI),
            lduAddr().patchSortStartAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}

namespace Foam
{
    template<class Type>
    struct fvMatrixAddBoundarySourceFunctor
    {
        const Type* pbc;
        const Type* pnf;
        const label* neiStart;
        const label* losort;

        fvMatrixAddBoundarySourceFunctor
        (
            const Type* _pbc,
            const Type* _pnf,
            const label* _neiStart,
            const label* _losort
        ):
             pbc(_pbc),
             pnf(_pnf),
             neiStart(_neiStart),
             losort(_losort)
        {}

        __host__ __device__
        Type operator()(const Type& d, const label& id)
        {
            Type out = d;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];

                out += cmptMultiply(pbc[face], pnf[face]);
            }

            return out;
        }
    };
}


template<class Type>
void Foam::fvMatrix<Type>::addBoundarySource
(
    gpuField<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const gpuField<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField
            (
                lduAddr().patchSortCells(patchI),
                lduAddr().patchSortAddr(patchI),
                lduAddr().patchSortStartAddr(patchI),
                pbc, 
                source
            );
        }
        else if (couples)
        {
            tmp<gpuField<Type> > tpnf = ptf.patchNeighbourField();
            const gpuField<Type>& pnf = tpnf();

            const labelgpuList& addr = lduAddr().patchSortCells(patchI);
            const labelgpuList& sort = lduAddr().patchSortAddr(patchI);
            const labelgpuList& sortStart = lduAddr().patchSortStartAddr(patchI);

            thrust::transform
            (
                thrust::make_permutation_iterator
                (
                    source.begin(),
                    addr.begin()
                ),
                thrust::make_permutation_iterator
                (
                    source.begin(),
                    addr.end()
                ),
                thrust::make_counting_iterator(0),
                thrust::make_permutation_iterator
                (
                    source.begin(),
                    addr.begin()
                ),
                fvMatrixAddBoundarySourceFunctor<Type>
                (
                    pbc.data(),
                    pnf.data(),
                    sortStart.data(),
                    sort.data()
                )
            );
        }
    }
}


namespace Foam
{
template<class Type>
struct fvMatrixSetValuesSourceFunctor : public std::binary_function<scalar,label,scalar>
{
    const bool* ownMask;
    const bool* neiMask;
    const Type* value;
    const scalar* upper;
    const scalar* lower;
    const label* ownStart;
    const label* neiStart;
    const label* own;
    const label* nei;
    const label* losort;

    fvMatrixSetValuesSourceFunctor
    (
        const bool* _ownMask,
        const bool* _neiMask,
        const Type* _value,
        const scalar* _upper,
        const scalar* _lower,
        const label* _ownStart,
        const label* _neiStart,
        const label* _own,
        const label* _nei,
        const label* _losort
    ):
        ownMask(_ownMask),
        neiMask(_neiMask),
        value(_value),
        upper(_upper),
        lower(_lower),
        ownStart(_ownStart),
        neiStart(_neiStart),
        own(_own),
        nei(_nei),
        losort(_losort)
    {}

    __HOST____DEVICE__
    Type operator()(const Type& source,const thrust::tuple<label,bool>& t)
    {
        Type out = source;
        label id = thrust::get<0>(t);
        bool cellSet = thrust::get<1>(t);
		
        if( ! cellSet)
        {
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;

            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                bool neiSet = neiMask[face];
				
                if(neiSet)
                {
                    out -= lower[face]*value[nei[face]];
                }
            }

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                bool ownSet = ownMask[face];
				
                if(ownSet)
                {
                    out -= upper[face]*value[own[face]];
                }
            }
        }

        return out;
    }
};

template<class Type>
struct fvMatrixSetValuesClearFacesFunctor
{
    const Type zero;

    fvMatrixSetValuesClearFacesFunctor(Type _zero): zero(_zero) {}

    __HOST____DEVICE__
    Type operator()(const Type& val,const bool& set)
    {
        if(set)
            return zero;
        else
            return val;
    }
};

}

template<class Type>
template<template<class> class ListType>
void Foam::fvMatrix<Type>::setValuesFromList
(
    const labelgpuList& cellLabels,
    const ListType<Type>& values
)
{
    const fvMesh& mesh = psi_.mesh();

    const labelgpuList& own = mesh.owner();
    const labelgpuList& nei = mesh.neighbour();

    scalargpuField& Diag = diag();
    gpuField<Type>& psi =
        const_cast
        <
            GeometricField<Type, fvPatchField, volMesh>&
        >(psi_).internalField();

    thrust::copy
    (
        values.begin(),
        values.end(),
        thrust::make_permutation_iterator
        (
            psi.begin(),
            cellLabels.begin()
        )
    );

    thrust::transform
    (
        values.begin(),
        values.end(),
        thrust::make_permutation_iterator
        (
            Diag.begin(),
            cellLabels.begin()
        ),
        thrust::make_permutation_iterator
        (
            source_.begin(),
            cellLabels.begin()
        ),
        multiplyOperatorFunctor<Type,scalar,Type>()
    );

    if (symmetric() || asymmetric())
    {
        //TODO make it somehow more comprehensible
        gpuList<bool> cellMask(Diag.size(),false);
        gpuList<bool> ownMask(own.size(),false);
        gpuList<bool> neiMask(nei.size(),false);
		
        thrust::fill
        (
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                cellLabels.begin()
            ),
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                cellLabels.end()
            ),
            true
        );
		             
        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                own.begin()
            ),
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                own.end()
            ),
            ownMask.begin()
        );
		             
        thrust::copy
        (
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                nei.begin()
            ),
            thrust::make_permutation_iterator
            (
                cellMask.begin(),
                nei.end()
            ),
            neiMask.begin()
        );
		             
        gpuList<Type> cellValuesTmp(Diag.size(),pTraits<Type>::zero);
		
        thrust::copy
        (
            values.begin(),
            values.end(),
            thrust::make_permutation_iterator
            (
                cellValuesTmp.begin(),
                cellLabels.begin()
            )
        );
		             
        const labelgpuList& l = lduAddr().lowerAddr();
        const labelgpuList& u = lduAddr().upperAddr();
        const labelgpuList& losort = lduAddr().losortAddr();

        const labelgpuList& ownStart = lduAddr().ownerStartAddr();
        const labelgpuList& losortStart = lduAddr().losortStartAddr();

        const scalargpuField& Lower = lower();
        const scalargpuField& Upper = upper();
		
        thrust::transform
        (
            source_.begin(),
            source_.end(),
            thrust::make_zip_iterator(thrust::make_tuple
            (
                thrust::make_counting_iterator(0),
                cellMask.begin()
            )),
            source_.begin(),
            fvMatrixSetValuesSourceFunctor<Type>
            (
                ownMask.data(),
                neiMask.data(),
                cellValuesTmp.data(),
                Upper.data(),
                Lower.data(),
                ownStart.data(),
                losortStart.data(),
                l.data(),
                u.data(),
                losort.data()
            )
        );


        thrust::transform
        (
            upper().begin(),
            upper().end(),
            ownMask.begin(),
            upper().begin(),
            fvMatrixSetValuesClearFacesFunctor<scalar>(0.0)
        );
                          
        if (asymmetric())
        {
            thrust::transform
            (
                lower().begin(),
                lower().end(),
                neiMask.begin(),
                lower().begin(),
                fvMatrixSetValuesClearFacesFunctor<scalar>(0.0)
            );
        }
        
        forAll(mesh.boundaryMesh(),patchi)
        {
            const gpuField<Type>& internalCoeffs = internalCoeffs_[patchi];
            const gpuField<Type>& boundaryCoeffs = boundaryCoeffs_[patchi];
            const labelgpuList& pcells = mesh.boundary()[patchi].faceCells();
			
            thrust::transform
            (
                internalCoeffs.begin(),
                internalCoeffs.end(),
                thrust::make_permutation_iterator
                (
                    cellMask.begin(),
                    pcells.begin()
                ),
                const_cast<gpuField<Type>&>(internalCoeffs).begin(),
                fvMatrixSetValuesClearFacesFunctor<Type>(pTraits<Type>::zero)
            );
                              
            thrust::transform
            (
                boundaryCoeffs.begin(),
                boundaryCoeffs.end(),
                thrust::make_permutation_iterator
                (
                    cellMask.begin(),
                    pcells.begin()
                ),
                const_cast<gpuField<Type>&>(boundaryCoeffs).begin(),
                fvMatrixSetValuesClearFacesFunctor<Type>(pTraits<Type>::zero)
            );
        }
        
    }
/*
    forAll(cellLabels, i)
    {
        const label celli = cellLabels[i];
        const Type& value = values[i];

//        psi[celli] = value;
//        source_[celli] = value*Diag[celli];

        if (symmetric() || asymmetric())
        {
            const cell& c = cells[celli];

            forAll(c, j)
            {
                const label facei = c[j];

                if (mesh.isInternalFace(facei))
                {
                    if (symmetric())
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= upper()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                    }
                    else
                    {
                        if (celli == own[facei])
                        {
                            source_[nei[facei]] -= lower()[facei]*value;
                        }
                        else
                        {
                            source_[own[facei]] -= upper()[facei]*value;
                        }

                        upper()[facei] = 0.0;
                        lower()[facei] = 0.0;
                    }
                }
                else
                {
                    label patchi = mesh.boundaryMesh().whichPatch(facei);

                    if (internalCoeffs_[patchi].size())
                    {
                        label patchFacei =
                            mesh.boundaryMesh()[patchi].whichFace(facei);

                        internalCoeffs_[patchi][patchFacei] =
                            pTraits<Type>::zero;

                        boundaryCoeffs_[patchi][patchFacei] =
                            pTraits<Type>::zero;
                    }
                }
            }
        }
    }
*/
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), pTraits<Type>::zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>(GeometricField<Type, fvPatchField, volMesh>&,"
               " const dimensionSet&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new gpuField<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new gpuField<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    GeometricField<Type, fvPatchField, volMesh>& psiRef =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryField().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const fvMatrix<Type>& fvm)
:
    refCount(),
    lduMatrix(fvm),
    psi_(fvm.psi_),
    dimensions_(fvm.dimensions_),
    source_(fvm.source_),
    internalCoeffs_(fvm.internalCoeffs_),
    boundaryCoeffs_(fvm.boundaryCoeffs_),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::fvMatrix(const fvMatrix<Type>&) : "
            << "copying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (fvm.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *(fvm.faceFluxCorrectionPtr_)
        );
    }
}


#ifdef ConstructFromTmp
template<class Type>
Foam::fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type> >& tfvm)
:
    refCount(),
    lduMatrix
    (
        const_cast<fvMatrix<Type>&>(tfvm()),
        tfvm.isTmp()
    ),
    psi_(tfvm().psi_),
    dimensions_(tfvm().dimensions_),
    source_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).source_,
        tfvm.isTmp()
    ),
    internalCoeffs_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).internalCoeffs_,
        tfvm.isTmp()
    ),
    boundaryCoeffs_
    (
        const_cast<fvMatrix<Type>&>(tfvm()).boundaryCoeffs_,
        tfvm.isTmp()
    ),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::fvMatrix(const tmp<fvMatrix<Type> >&) : "
            << "copying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (tfvm().faceFluxCorrectionPtr_)
    {
        if (tfvm.isTmp())
        {
            faceFluxCorrectionPtr_ = tfvm().faceFluxCorrectionPtr_;
            tfvm().faceFluxCorrectionPtr_ = NULL;
        }
        else
        {
            faceFluxCorrectionPtr_ = new
                GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    *(tfvm().faceFluxCorrectionPtr_)
                );
        }
    }

    tfvm.clear();
}
#endif


template<class Type>
Foam::fvMatrix<Type>::fvMatrix
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>"
               "(GeometricField<Type, fvPatchField, volMesh>&, Istream&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new gpuField<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new gpuField<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

}


template<class Type>
Foam::fvMatrix<Type>::~fvMatrix()
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::~fvMatrix<Type>() : "
            << "destroying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelgpuList& cellLabels,
    const gpuList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}

/*
template<class Type>
void Foam::fvMatrix<Type>::setValues
(
    const labelgpuList& cellLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}
*/
template<class Type>
void Foam::fvMatrix<Type>::setReference
(
    const label celli,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && celli >= 0)
    {
        source().set(celli,source().get(celli)+diag().get(celli)*value);
        diag().set(celli,2*diag().get(celli));
/*
        source()[celli] += diag()[celli]*value;
        diag()[celli] += diag()[celli];
*/
    }
}

namespace Foam
{
    template<class Type>
    struct fvMatrixRelaxDiagonalDominanceFunctor
    {
        __HOST____DEVICE__
        scalar operator()(const scalar& s1, const scalar& s2)
        {
            return max(mag(s1), s2);
        }
    };


    template<class Type,class Fun>
    struct fvMatrixRelaxAddToDiagonalFunctor : public std::binary_function<label,scalar,scalar>
    {
        const Fun f;
        const Type* iCoeffs;
        const label* neiStart;
        const label* losort;

        fvMatrixRelaxAddToDiagonalFunctor
        (
            const Fun _f,
            const Type* _iCoeffs,
            const label* _neiStart,
            const label* _losort
        ):
             f(_f),
             iCoeffs(_iCoeffs),
             neiStart(_neiStart),
             losort(_losort)
        {}

        __HOST____DEVICE__
        scalar operator()(const label& id,const scalar& s)
        {
            scalar out = s;

            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;

            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out += f(iCoeffs[face]);
            }

            return out;
        }
    };
    
    template<class Type>
    struct componetZeroFunctor
    {
        __HOST____DEVICE__
        scalar operator ()(const Type& t) const
        {
            return component(t, 0);
        }
    };
	
    template<class Type>
    struct negativeComponetZeroFunctor
    {
        __HOST____DEVICE__
        scalar operator ()(const Type& t) const
        {
            return -component(t, 0);
        }
    };
	
    template<class Type>
    struct magComponetZeroFunctor
    {
        __HOST____DEVICE__
        scalar operator ()(const Type& t) const
        {
            return mag(component(t, 0));
        }
    };
	
    template<class Type>
    struct maxComponentMagComponetFunctor
    {
        __HOST____DEVICE__
        scalar operator ()(const Type& t) const
        {
            return cmptMax(cmptMag(t));
        }
    };
	
    template<class Type>
    struct negativeComponetMinFunctor
    {
        __HOST____DEVICE__
        scalar operator ()(const Type& t) const
        {
            return - cmptMin(t);
        }
    };
}


template<class Type>
void Foam::fvMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    if (debug)
    {
        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
            << "Relaxing " << psi_.name() << " by " << alpha
            << endl;
    }

    gpuField<Type>& S = source();
    scalargpuField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalargpuField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalargpuField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            gpuField<Type>& iCoeffs = internalCoeffs_[patchI];
            
            const labelgpuList& pcells = lduAddr().patchSortCells(patchI);

            const labelgpuList& losort = lduAddr().patchSortAddr(patchI);
            const labelgpuList& losortStart = lduAddr().patchSortStartAddr(patchI);

            if (ptf.coupled())
            {
                const gpuField<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                                  
                thrust::transform
                (
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0)+pcells.size(),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    fvMatrixRelaxAddToDiagonalFunctor<Type,componetZeroFunctor<Type> >
                    (
                        componetZeroFunctor<Type>(),
                        iCoeffs.data(),
                        losortStart.data(),
                        losort.data()
                    )
                );
                                               
                thrust::transform
                (
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0)+pcells.size(),
                    thrust::make_permutation_iterator
                    (
                        sumOff.begin(),
                        pcells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        sumOff.begin(),
                        pcells.begin()
                    ),
                    fvMatrixRelaxAddToDiagonalFunctor<Type,magComponetZeroFunctor<Type> >
                    (
                        magComponetZeroFunctor<Type>(),
                        pCoeffs.data(),
                        losortStart.data(),
                        losort.data()
                    )
                );

            }
            else
            {
                // For non-coupled boundaries add the maximum magnitude diagonal
                // contribution to ensure stability
                                  
                thrust::transform
                (
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0)+pcells.size(),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    fvMatrixRelaxAddToDiagonalFunctor<Type,maxComponentMagComponetFunctor<Type> >
                    (
                        maxComponentMagComponetFunctor<Type>(),
                        iCoeffs.data(),
                        losortStart.data(),
                        losort.data()
                    )
                );
            }
        }
    }

/* 
    if (debug)
    {
        // Calculate amount of non-dominance.
        label nNon = 0;
        scalar maxNon = 0.0;
        scalar sumNon = 0.0;
        forAll(D, celli)
        {
            scalar d = (sumOff[celli] - D[celli])/mag(D[celli]);

            if (d > 0)
            {
                nNon++;
                maxNon = max(maxNon, d);
                sumNon += d;
            }
        }

        reduce(nNon, sumOp<label>(), UPstream::msgType(), psi_.mesh().comm());
        reduce
        (
            maxNon,
            maxOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        reduce
        (
            sumNon,
            sumOp<scalar>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );
        sumNon /= returnReduce
        (
            D.size(),
            sumOp<label>(),
            UPstream::msgType(),
            psi_.mesh().comm()
        );

        InfoIn("fvMatrix<Type>::relax(const scalar alpha)")
            << "Matrix dominance test for " << psi_.name() << nl
            << "    number of non-dominant cells   : " << nNon << nl
            << "    maximum relative non-dominance : " << maxNon << nl
            << "    average relative non-dominance : " << sumNon << nl
            << endl;
    }
*/

    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
/*    forAll(D, celli)
    {
        D[celli] = max(mag(D[celli]), sumOff[celli]);
    }*/

    thrust::transform(D.begin(),D.end(),sumOff.begin(),D.begin(),
                      fvMatrixRelaxDiagonalDominanceFunctor<Type>());

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            gpuField<Type>& iCoeffs = internalCoeffs_[patchI];
            
            const labelgpuList& pcells = lduAddr().patchSortCells(patchI);

            const labelgpuList& losort = lduAddr().patchSortAddr(patchI);
            const labelgpuList& losortStart = lduAddr().patchSortStartAddr(patchI);

            if (ptf.coupled())
            {
                thrust::transform
                (
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0)+pcells.size(),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    fvMatrixRelaxAddToDiagonalFunctor<Type,negativeComponetZeroFunctor<Type> >(
                        negativeComponetZeroFunctor<Type>(),
                        iCoeffs.data(),
                        losortStart.data(),
                        losort.data()
                    )
                );
            }
            else
            {              
                thrust::transform
                (
                    thrust::make_counting_iterator(0),
                    thrust::make_counting_iterator(0)+pcells.size(),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    thrust::make_permutation_iterator
                    (
                        D.begin(),
                        pcells.begin()
                    ),
                    fvMatrixRelaxAddToDiagonalFunctor<Type,negativeComponetMinFunctor<Type> >
                    (
                        negativeComponetMinFunctor<Type>(),
                        iCoeffs.data(),
                        losortStart.data(),
                        losort.data()
                    )
                );

            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.internalField();
}


template<class Type>
void Foam::fvMatrix<Type>::relax()
{
    word name = psi_.select
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
    );

    if (psi_.mesh().relaxEquation(name))
    {
        relax(psi_.mesh().equationRelaxationFactor(name));
    }
}


template<class Type>
void Foam::fvMatrix<Type>::boundaryManipulate
(
    typename GeometricField<Type, fvPatchField, volMesh>::
        GeometricBoundaryField& bFields
)
{
    forAll(bFields, patchI)
    {
        bFields[patchI].manipulateMatrix(*this);
    }
}

template<class Type>
void Foam::fvMatrix<Type>::D(Foam::scalargpuField& tdiag) const
{
    tdiag = diag();
    addCmptAvBoundaryDiag(tdiag);
}

template<class Type>
Foam::tmp<Foam::scalargpuField> Foam::fvMatrix<Type>::D() const
{
    tmp<scalargpuField> tdiag(new scalargpuField(diag().size()));
    D(tdiag());
    return tdiag;
}

template<class Type>
void Foam::fvMatrix<Type>::DD(Foam::gpuField<Type>& tdiag) const
{
    tdiag = pTraits<Type>::one*diag();

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchSortCells(patchI),
                lduAddr().patchSortAddr(patchI),
                lduAddr().patchSortStartAddr(patchI),
                internalCoeffs_[patchI],
                tdiag
            );
        }
    }
}

template<class Type>
Foam::tmp<Foam::gpuField<Type> > Foam::fvMatrix<Type>::DD() const
{
    tmp<gpuField<Type> > tdiag(diag().size());

    DD(tdiag());

    return tdiag;
}

template<class Type>
void Foam::fvMatrix<Type>::A(Foam::volScalarField& Aphi) const
{
    Aphi.internalField() = D()/psi_.mesh().V().getField();
    Aphi.correctBoundaryConditions();
}

template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::A() const
{
    tmp<volScalarField> tAphi
    (
        new volScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    A(tAphi());

    return tAphi;
}

template<class Type>
void Foam::fvMatrix<Type>::H(Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& Hphi) const
{
    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        label pSize = psi_.size();

        scalargpuField psiCmpt(fvMatrixCache::first(level(),pSize),pSize);
        component(psiCmpt,psi_.internalField(),cmpt);

        scalargpuField boundaryDiagCmpt(fvMatrixCache::second(level(),pSize),pSize);
        boundaryDiagCmpt = 0.0;

        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.internalField().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    lduMatrix::H(Hphi.internalField(),psi_.internalField());

    Hphi.internalField() += source_;
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().V().getField();
    Hphi.correctBoundaryConditions();

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1)
        {
            Hphi.replace
            (
                cmpt,
                dimensionedScalar("0", Hphi.dimensions(), 0.0)
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::fvMatrix<Type>::H() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tHphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Hphi = tHphi();

    H(Hphi);

    return tHphi;
}

template<class Type>
void Foam::fvMatrix<Type>::H1(Foam::volScalarField& H1_) const
{
    lduMatrix::H1(H1_.internalField());

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchSortCells(patchI),
                lduAddr().patchSortAddr(patchI),
                lduAddr().patchSortStartAddr(patchI),
                boundaryCoeffs_[patchI].component(0),
                H1_.getField()
            );
        }
    }

    H1_.internalField() /= psi_.mesh().V().getField();
    H1_.correctBoundaryConditions();
}

template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Type>::H1() const
{
    tmp<volScalarField> tH1
    (
        new volScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1();

    H1(H1_);

    return tH1;
}

template<class Type>
void Foam::fvMatrix<Type>::flux
(
    Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>& fieldFlux
) const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorIn("fvMatrix<Type>::flux()")
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        label pSize = psi_.size();

        scalargpuField faceHTmp(fvMatrixCache::first(level(),lower().size()),lower().size());
        scalargpuField psiTmp(fvMatrixCache::second(level(),pSize),pSize);

        component(psiTmp,psi_.internalField(),cmpt);
        lduMatrix::faceH(faceHTmp,psiTmp);

        fieldFlux.internalField().replace
        (
            cmpt,
            faceHTmp
        );
    }

    FieldField<gpuField, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<gpuField, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryField()[patchI] =
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fvMatrix<Type>::
flux() const
{
    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfieldFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& fieldFlux = tfieldFlux();

    flux(fieldFlux);

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::operator=(const fvMatrix<Type>& fvmv)
{
    if (this == &fvmv)
    {
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
            << "different fields"
            << abort(FatalError);
    }

    dimensions_ = fvmv.dimensions_;
    lduMatrix::operator=(fvmv);
    source_ = fvmv.source_;
    internalCoeffs_ = fvmv.internalCoeffs_;
    boundaryCoeffs_ = fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
        (*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ += fvmv.dimensions_;
    lduMatrix::operator+=(fvmv);
    source_ += fvmv.source_;
    internalCoeffs_ += fvmv.internalCoeffs_;
    boundaryCoeffs_ += fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ -= fvmv.dimensions_;
    lduMatrix::operator-=(fvmv);
    source_ -= fvmv.source_;
    internalCoeffs_ -= fvmv.internalCoeffs_;
    boundaryCoeffs_ -= fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
        (-*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator-=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().V().getField()*su.getField();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().V().getField()*su.getField();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= psi().mesh().V().getField()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += psi().mesh().V().getField()*su;
}


template<class Type>
void Foam::fvMatrix<Type>::operator+=
(
    const zero&
)
{}


template<class Type>
void Foam::fvMatrix<Type>::operator-=
(
    const zero&
)
{}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const DimensionedField<scalar, volMesh>& dsf
)
{
    dimensions_ *= dsf.dimensions();
    lduMatrix::operator*=(dsf.getField());
    source_ *= dsf.getField();

    forAll(boundaryCoeffs_, patchI)
    {
        scalargpuField psf
        (
            dsf.mesh().boundary()[patchI].patchInternalField(dsf.getField())
        );

        internalCoeffs_[patchI] *= psf;
        boundaryCoeffs_[patchI] *= psf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::operator*="
            "(const DimensionedField<scalar, volMesh>&)"
        )   << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<DimensionedField<scalar, volMesh> >& tdsf
)
{
    operator*=(tdsf());
    tdsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const volScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.getField());
    source_ *= vsf.getField();

    forAll(vsf.boundaryField(), patchI)
    {
        const fvPatchScalarField& psf = vsf.boundaryField()[patchI];

        if (psf.coupled())
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf.patchNeighbourField();
        }
        else
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf;
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::operator*="
            "(const DimensionedField<scalar, volMesh>&)"
        )   << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}

template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const tmp<volScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void Foam::fvMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm1,
    const fvMatrix<Type>& fvm2,
    const char* op
)
{
    if (&fvm1.psi() != &fvm2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << fvm1.dimensions()/dimVolume << " ] "
            << op
            << " [" << fvm2.psi().name() << fvm2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const DimensionedField<Type, volMesh>& df,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != df.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const GeometricField<Type, "
            "fvPatchField, volMesh>&)"
        )   <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << df.name() << df.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const fvMatrix<Type>& fvm,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::solverPerformance Foam::solve
(
    fvMatrix<Type>& fvm,
    const dictionary& solverControls
)
{
    return fvm.solve(solverControls);
}

template<class Type>
Foam::solverPerformance Foam::solve
(
    const tmp<fvMatrix<Type> >& tfvm,
    const dictionary& solverControls
)
{
    solverPerformance solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve(solverControls);

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::solverPerformance Foam::solve(fvMatrix<Type>& fvm)
{
    return fvm.solve();
}

template<class Type>
Foam::solverPerformance Foam::solve(const tmp<fvMatrix<Type> >& tfvm)
{
    solverPerformance solverPerf =
        const_cast<fvMatrix<Type>&>(tfvm()).solve();

    tfvm.clear();

    return solverPerf;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::correction
(
    const fvMatrix<Type>& A
)
{
    tmp<Foam::fvMatrix<Type> > tAcorr = A - (A & A.psi());

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().fluxRequired(A.psi().name())
    )
    {
        tAcorr().faceFluxCorrectionPtr() = (-A.flux()).ptr();
    }

    return tAcorr;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::correction
(
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<Foam::fvMatrix<Type> > tAcorr = tA - (tA() & tA().psi());

    // Note the matrix coefficients are still that of matrix A
    const fvMatrix<Type>& A = tAcorr();

    if
    (
        (A.hasUpper() || A.hasLower())
     && A.psi().mesh().fluxRequired(A.psi().name())
    )
    {
        tAcorr().faceFluxCorrectionPtr() = (-A.flux()).ptr();
    }

    return tAcorr;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += A.psi().mesh().V().getField()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tC().psi().mesh().V().getField()*su.value();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const fvMatrix<Type>& A,
    const zero&
)
{
    return A;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const zero&
)
{
    return tA;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const DimensionedField<Type, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<DimensionedField<Type, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const DimensionedField<Type, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.mesh().V().getField()*su.getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V().getField()*tsu().getField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V().getField()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.value()*tC().psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.value()*tC().psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.value()*A.psi().mesh().V().getField();
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.value()*tC().psi().mesh().V().getField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const DimensionedField<scalar, volMesh>& dsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp< DimensionedField<scalar, volMesh> >& tdsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const DimensionedField<scalar, volMesh>& dsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= dsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<DimensionedField<scalar, volMesh> >& tdsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= tdsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= ds;
    return tC;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> > Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "M&" + psi.name(),
                psi.instance(),
                psi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi.mesh(),
            M.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Mphi = tMphi();

    // Loop over field components
    if (M.hasDiag())
    {
        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            scalargpuField psiCmpt(psi.getField().component(cmpt));
            scalargpuField boundaryDiagCmpt(M.diag());
            M.addBoundaryDiag(boundaryDiagCmpt, cmpt);
            Mphi.internalField().replace(cmpt, -boundaryDiagCmpt*psiCmpt);
        }
    }
    else
    {
        Mphi.internalField() = pTraits<Type>::zero;
    }

    Mphi.internalField() += M.lduMatrix::H(psi.getField()) + M.source();
    M.addBoundarySource(Mphi.internalField());

    Mphi.internalField() /= -psi.mesh().V().getField();
    Mphi.correctBoundaryConditions();

    return tMphi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const fvMatrix<Type>& M,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::operator&
(
    const tmp<fvMatrix<Type> >& tM,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvMatrix<Type>& fvm)
{
    os  << static_cast<const lduMatrix&>(fvm) << nl
        << fvm.dimensions_ << nl
        << fvm.source_ << nl
        << fvm.internalCoeffs_ << nl
        << fvm.boundaryCoeffs_ << endl;

    os.check("Ostream& operator<<(Ostream&, fvMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "fvMatrixSolve.C"

// ************************************************************************* //
