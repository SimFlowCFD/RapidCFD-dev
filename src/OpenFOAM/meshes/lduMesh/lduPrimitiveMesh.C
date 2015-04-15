/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "lduPrimitiveMesh.H"
#include "processorLduInterface.H"
#include "EdgeMap.H"
#include "labelPair.H"
#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMesh, 0);

    //- Less operator for pairs of \<processor\>\<index\>
    class procLess
    {
        const labelPairList& lst_;

    public:

        procLess(const labelPairList& lst)
        :
            lst_(lst)
        {}

        bool operator()(const label a, const label b)
        {
            return lst_[a].first() < lst_[b].first();
        }
    };
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduPrimitiveMesh::checkUpperTriangular
(
    const label size,
    const labelUList& l,
    const labelUList& u
)
{
    forAll(l, faceI)
    {
        if (u[faceI] < l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Reversed face. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << abort(FatalError);
        }
        if (l[faceI] < 0 || u[faceI] < 0 || u[faceI] >= size)
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Illegal cell label. Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << abort(FatalError);
        }
    }

    for (label faceI=1; faceI < l.size(); faceI++)
    {
        if (l[faceI-1] > l[faceI])
        {
            FatalErrorIn
            (
                "checkUpperTriangular"
                "(const label, const labelUList&, const labelUList&)"
            )   << "Lower not in incremental cell order."
                << " Problem at face " << faceI
                << " l:" << l[faceI] << " u:" << u[faceI]
                << " previous l:" << l[faceI-1]
                << abort(FatalError);
        }
        else if (l[faceI-1] == l[faceI])
        {
            // Same cell.
            if (u[faceI-1] > u[faceI])
            {
                FatalErrorIn
                (
                    "checkUpperTriangular"
                    "(const label, const labelUList&, const labelUList&)"
                )   << "Upper not in incremental cell order."
                    << " Problem at face " << faceI
                    << " l:" << l[faceI] << " u:" << u[faceI]
                    << " previous u:" << u[faceI-1]
                    << abort(FatalError);
            }
        }
    }
}


Foam::label Foam::lduPrimitiveMesh::totalSize
(
    const PtrList<lduPrimitiveMesh>& meshes
)
{
    label size = 0;

    forAll(meshes, i)
    {
        size += meshes[i].lduAddr().size();
    }
    return size;
}

Foam::labelList Foam::lduPrimitiveMesh::upperTriOrder
(
    const label nCells,
    const labelUList& lower,
    const labelUList& upper
)
{
    labelList nNbrs(nCells, 0);

    // Count number of upper neighbours
    forAll(lower, faceI)
    {
        if (upper[faceI] < lower[faceI])
        {
            FatalErrorIn("lduPrimitiveMesh::upperTriOrder(..)")
                << "Problem at face:" << faceI
                << " lower:" << lower[faceI]
                << " upper:" << upper[faceI]
                << exit(FatalError);
        }
        nNbrs[lower[faceI]]++;
    }

    // Construct cell-upper cell addressing
    labelList offsets(nCells+1);
    offsets[0] = 0;
    forAll(nNbrs, cellI)
    {
        offsets[cellI+1] = offsets[cellI]+nNbrs[cellI];
    }

    nNbrs = offsets;

    labelList cellToFaces(offsets.last());
    forAll(upper, faceI)
    {
        label cellI = lower[faceI];
        cellToFaces[nNbrs[cellI]++] = faceI;
    }

    // Sort

    labelList oldToNew(lower.size());

    labelList order;
    labelList nbr;

    label newFaceI = 0;

    for (label cellI = 0; cellI < nCells; cellI++)
    {
        label startOfCell = offsets[cellI];
        label nNbr = offsets[cellI+1] - startOfCell;

        nbr.setSize(nNbr);
        order.setSize(nNbr);
        forAll(order, i)
        {
            nbr[i] = upper[cellToFaces[offsets[cellI]+i]];
        }
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            oldToNew[cellToFaces[startOfCell + index]] = newFaceI++;
        }
    }

    return oldToNew;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label level,
    const label nCells,
    labelgpuList& l,
    labelgpuList& u,
    const label comm,
    bool reUse
)
:
    lduAddressing(nCells),
    level_(level),
    lowerAddrHost_(l.size()),
    upperAddrHost_(u.size()),
    lowerAddr_(l, reUse),
    upperAddr_(u, reUse),
    comm_(comm)
{
    thrust::copy(lowerAddr_.begin(),lowerAddr_.end(),lowerAddrHost_.begin());
    thrust::copy(upperAddr_.begin(),upperAddr_.end(),upperAddrHost_.begin());
}

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label level,
    const label nCells,
    labelList& l,
    labelList& u,
    const label comm,
    bool reUse
)
:
    lduAddressing(nCells),
    level_(level),
    lowerAddrHost_(l, reUse),
    upperAddrHost_(u, reUse),
    comm_(comm)
{
    lowerAddr_=lowerAddrHost_;
    upperAddr_=upperAddrHost_;
}


void Foam::lduPrimitiveMesh::addInterfaces
(
    lduInterfacePtrsList& interfaces,
    const lduSchedule& ps
)
{
    interfaces_ = interfaces;

    // Create interfaces
    primitiveInterfaces_.setSize(interfaces_.size());
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            primitiveInterfaces_.set(i, &interfaces_[i]);
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label level,
    const label nCells,
    labelgpuList& l,
    labelgpuList& u,
    PtrList<const lduInterface>& primitiveInterfaces,
    const lduSchedule& ps,
    const label comm
)
:
    lduAddressing(nCells), 
    level_(level),   
    lowerAddrHost_(l.size()),
    upperAddrHost_(u.size()),
    lowerAddr_(l, true),
    upperAddr_(u, true),
    primitiveInterfaces_(0),
    patchSchedule_(ps),
    comm_(comm)
{
    thrust::copy(lowerAddr_.begin(),lowerAddr_.end(),lowerAddrHost_.begin());
    thrust::copy(upperAddr_.begin(),upperAddr_.end(),upperAddrHost_.begin());

    primitiveInterfaces_.transfer(primitiveInterfaces);

    // Create interfaces
    interfaces_.setSize(primitiveInterfaces_.size());
    forAll(primitiveInterfaces_, i)
    {
        if (primitiveInterfaces_.set(i))
        {
            interfaces_.set(i, &primitiveInterfaces_[i]);
        }
    }
}

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label level,
    const label nCells,
    labelList& l,
    labelList& u,
    PtrList<const lduInterface>& primitiveInterfaces,
    const lduSchedule& ps,
    const label comm
)
:
    lduAddressing(nCells),  
    level_(level),  
    lowerAddrHost_(l, true),
    upperAddrHost_(u, true),
    lowerAddr_(lowerAddrHost_),
    upperAddr_(upperAddrHost_),
    primitiveInterfaces_(0),
    patchSchedule_(ps),
    comm_(comm)
{
    primitiveInterfaces_.transfer(primitiveInterfaces);

    // Create interfaces
    interfaces_.setSize(primitiveInterfaces_.size());
    forAll(primitiveInterfaces_, i)
    {
        if (primitiveInterfaces_.set(i))
        {
            interfaces_.set(i, &primitiveInterfaces_[i]);
        }
    }
}

// ************************************************************************* //
