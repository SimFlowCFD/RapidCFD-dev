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

#include "GAMGSolver.H"
#include "GAMGInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "processorGAMGInterfaceField.H"
#include "gpuIndirectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

struct luGAMGNonNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x>=0;
    }
};

struct luGAMGNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x<0;
    }
};

struct faceToDiagFunctor 
   : public std::unary_function<label,label>
{
    __HOST____DEVICE__
    label operator()(const label& n)
    {
        return -1 - n;
    }
};


struct GAMGSolverAgglomerateAsymFunctor
{
    const scalar* uf;
    const scalar* lf;
    const bool* flip;
    const label* sort;

    GAMGSolverAgglomerateAsymFunctor
    (
        const scalar* _uf,
        const scalar* _lf,
        const bool* _flip,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        flip(_flip),
        sort(_sort)
    {}
    
    template<class Tuple,class Tuple2>
    __HOST____DEVICE__
    thrust::tuple<scalar,scalar> operator()(const Tuple& in,const Tuple2& t)
    {
        scalar uc = thrust::get<0>(in);
        scalar lc = thrust::get<1>(in); 
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            bool flipVal = flip[index];
            scalar ufVal = uf[index];
            scalar lfVal = lf[index];
            if( ! flipVal)
            {
                uc += ufVal;
                lc += lfVal;
            }
            else
            {
                uc += lfVal;
                lc += ufVal;
            }
        }

        return thrust::make_tuple(uc,lc);
    }

};


struct GAMGSolverAgglomerateDiagAsymFunctor
{
    const scalar* uf;
    const scalar* lf;
    const label* sort;

    GAMGSolverAgglomerateDiagAsymFunctor
    (
        const scalar* _uf,
        const scalar* _lf,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        sort(_sort)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            out += uf[index] + lf[index];
        }

        return out;
    }
};


struct GAMGSolverAgglomerateSymFunctor
{
    const scalar* ff;
    const label* sort;

    GAMGSolverAgglomerateSymFunctor
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}
    
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += ff[sort[i]];
        }

        return out;
    }

};

struct GAMGSolverAgglomerateDiagSymFunctor
{
    const scalar* ff;
    const label* sort;

    GAMGSolverAgglomerateDiagSymFunctor
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += 2*ff[sort[i]];
        }

        return out;
    }
};


}

void Foam::GAMGSolver::agglomerateMatrix
(
    const label fineLevelIndex,
    const lduMesh& coarseMesh,
    const lduInterfacePtrsList& coarseMeshInterfaces
)
{
    // Get fine matrix
    const lduMatrix& fineMatrix = matrixLevel(fineLevelIndex);

    if (UPstream::myProcNo(fineMatrix.mesh().comm()) != -1)
    {
        const label nCoarseFaces = agglomeration_.nFaces(fineLevelIndex);
        const label nCoarseCells = agglomeration_.nCells(fineLevelIndex);

        // Set the coarse level matrix
        matrixLevels_.set
        (
            fineLevelIndex,
            new lduMatrix(coarseMesh)
        );
        lduMatrix& coarseMatrix = matrixLevels_[fineLevelIndex];


        // Coarse matrix diagonal initialised by restricting the finer mesh
        // diagonal. Note that we size with the cached coarse nCells and not
        // the actual coarseMesh size since this might be dummy when processor
        // agglomerating.
        scalargpuField& coarseDiag = coarseMatrix.diag(nCoarseCells);

        agglomeration_.restrictField
        (
            coarseDiag,
            fineMatrix.diag(),
            fineLevelIndex,
            false               // no processor agglomeration
        );

        // Get reference to fine-level interfaces
        const lduInterfaceFieldPtrsList& fineInterfaces =
            interfaceLevel(fineLevelIndex);

        // Create coarse-level interfaces
        primitiveInterfaceLevels_.set
        (
            fineLevelIndex,
            new PtrList<lduInterfaceField>(fineInterfaces.size())
        );

        PtrList<lduInterfaceField>& coarsePrimInterfaces =
            primitiveInterfaceLevels_[fineLevelIndex];

        interfaceLevels_.set
        (
            fineLevelIndex,
            new lduInterfaceFieldPtrsList(fineInterfaces.size())
        );

        lduInterfaceFieldPtrsList& coarseInterfaces =
            interfaceLevels_[fineLevelIndex];

        // Set coarse-level boundary coefficients
        interfaceLevelsBouCoeffs_.set
        (
            fineLevelIndex,
            new FieldField<gpuField, scalar>(fineInterfaces.size())
        );
        FieldField<gpuField, scalar>& coarseInterfaceBouCoeffs =
            interfaceLevelsBouCoeffs_[fineLevelIndex];

        // Set coarse-level internal coefficients
        interfaceLevelsIntCoeffs_.set
        (
            fineLevelIndex,
            new FieldField<gpuField, scalar>(fineInterfaces.size())
        );
        FieldField<gpuField, scalar>& coarseInterfaceIntCoeffs =
            interfaceLevelsIntCoeffs_[fineLevelIndex];

        // Add the coarse level
        agglomerateInterfaceCoefficients
        (
            fineLevelIndex,
            coarseMeshInterfaces,
            coarsePrimInterfaces,
            coarseInterfaces,
            coarseInterfaceBouCoeffs,
            coarseInterfaceIntCoeffs
        );


        // Get face restriction addressing for current level
        const labelgpuList& faceRestrictSortAddr =
            agglomeration_.faceRestrictSortAddressing(fineLevelIndex); 
        const labelgpuList& faceRestrictTargetAddr = 
            agglomeration_.faceRestrictTargetAddressing(fineLevelIndex);
        const labelgpuList& faceRestrictTargetStartAddr = 
            agglomeration_.faceRestrictTargetStartAddressing(fineLevelIndex);
        const boolgpuList& faceFlipMap =
            agglomeration_.faceFlipMap(fineLevelIndex);

        typedef typename thrust::pair<labelgpuList::iterator,scalargpuField::iterator> Pair;

        // Check if matrix is asymetric and if so agglomerate both upper
        // and lower coefficients ...
        if (fineMatrix.hasLower())
        {
            // Get off-diagonal matrix coefficients
            const scalargpuField& fineUpper = fineMatrix.upper();
            const scalargpuField& fineLower = fineMatrix.lower();

            // Coarse matrix upper coefficients. Note passed in size
            scalargpuField& coarseUpper = coarseMatrix.upper(nCoarseFaces);
            scalargpuField& coarseLower = coarseMatrix.lower(nCoarseFaces);


            thrust::transform_if
            (
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        coarseUpper.begin(),
                        faceRestrictTargetAddr.begin()
                    ),

                    thrust::make_permutation_iterator
                    (
                        coarseLower.begin(),
                        faceRestrictTargetAddr.begin()
                    )
                )),

                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        coarseUpper.begin(),
                        faceRestrictTargetAddr.end()
                    ),

                    thrust::make_permutation_iterator
                    (
                        coarseLower.begin(),
                        faceRestrictTargetAddr.end()
                    )
                )),

                thrust::make_zip_iterator(thrust::make_tuple
                (
                    faceRestrictTargetStartAddr.begin(),
                    faceRestrictTargetStartAddr.begin()+1
                )),

                faceRestrictTargetAddr.begin(),

                thrust::make_zip_iterator(thrust::make_tuple
                (
                    thrust::make_permutation_iterator
                    (
                        coarseUpper.begin(),
                        faceRestrictTargetAddr.begin()
                    ),

                    thrust::make_permutation_iterator
                    (
                        coarseLower.begin(),
                        faceRestrictTargetAddr.begin()
                    )
                )),

                GAMGSolverAgglomerateAsymFunctor
                (
                    fineUpper.data(),
                    fineLower.data(),
                    faceFlipMap.data(),
                    faceRestrictSortAddr.data()
                ),

                luGAMGNonNegative()
            );


            thrust::transform_if
            (
                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.begin(),
                        faceToDiagFunctor()
                    )
                ),

                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.end(),
                        faceToDiagFunctor()
                    )
                ),

                thrust::make_zip_iterator(thrust::make_tuple
                (
                    faceRestrictTargetStartAddr.begin(),
                    faceRestrictTargetStartAddr.begin()+1
                )),

                faceRestrictTargetAddr.begin(),

                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.begin(),
                        faceToDiagFunctor()
                    )
                ),

                GAMGSolverAgglomerateDiagAsymFunctor
                (
                    fineUpper.data(),
                    fineLower.data(),
                    faceRestrictSortAddr.data()
                ),

                luGAMGNegative()
            );

        }
        else // ... Otherwise it is symmetric so agglomerate just the upper
        {
            // Get off-diagonal matrix coefficients
            const scalargpuField& fineUpper = fineMatrix.upper();

            // Coarse matrix upper coefficients
            scalargpuField& coarseUpper = coarseMatrix.upper(nCoarseFaces);

            thrust::transform_if
            (
                thrust::make_permutation_iterator
                (
                    coarseUpper.begin(),
                    faceRestrictTargetAddr.begin()
                ),
                thrust::make_permutation_iterator
                (
                    coarseUpper.begin(),
                    faceRestrictTargetAddr.end()
                ),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    faceRestrictTargetStartAddr.begin(),
                    faceRestrictTargetStartAddr.begin()+1
                )),
                faceRestrictTargetAddr.begin(),
                thrust::make_permutation_iterator
                (
                    coarseUpper.begin(),
                    faceRestrictTargetAddr.begin()
                ),
                GAMGSolverAgglomerateSymFunctor
                (
                    fineUpper.data(),
                    faceRestrictSortAddr.data()
                ),
                luGAMGNonNegative()
            );

            thrust::transform_if
            (
                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.begin(),
                        faceToDiagFunctor()
                    )
                ),
                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.end(),
                        faceToDiagFunctor()
                    )
                ),
                thrust::make_zip_iterator(thrust::make_tuple
                (
                    faceRestrictTargetStartAddr.begin(),
                    faceRestrictTargetStartAddr.begin()+1
                )),
                faceRestrictTargetAddr.begin(),
                thrust::make_permutation_iterator
                (
                    coarseDiag.begin(),
                    thrust::make_transform_iterator
                    (
                        faceRestrictTargetAddr.begin(),
                        faceToDiagFunctor()
                    )
                ),
                GAMGSolverAgglomerateDiagSymFunctor
                (
                    fineUpper.data(),
                    faceRestrictSortAddr.data()
                ),
                luGAMGNegative()
            );
        }
    }
}


// Agglomerate only the interface coefficients.
void Foam::GAMGSolver::agglomerateInterfaceCoefficients
(
    const label fineLevelIndex,
    const lduInterfacePtrsList& coarseMeshInterfaces,
    PtrList<lduInterfaceField>& coarsePrimInterfaces,
    lduInterfaceFieldPtrsList& coarseInterfaces,
    FieldField<gpuField, scalar>& coarseInterfaceBouCoeffs,
    FieldField<gpuField, scalar>& coarseInterfaceIntCoeffs
) const
{
    // Get reference to fine-level interfaces
    const lduInterfaceFieldPtrsList& fineInterfaces =
        interfaceLevel(fineLevelIndex);

    // Get reference to fine-level boundary coefficients
    const FieldField<gpuField, scalar>& fineInterfaceBouCoeffs =
        interfaceBouCoeffsLevel(fineLevelIndex);

    // Get reference to fine-level internal coefficients
    const FieldField<gpuField, scalar>& fineInterfaceIntCoeffs =
        interfaceIntCoeffsLevel(fineLevelIndex);

    const labelgpuListList& patchFineToCoarseSort =
        agglomeration_.patchFaceRestrictSortAddressing(fineLevelIndex);

    const labelgpuListList& patchFineToCoarseTarget =
        agglomeration_.patchFaceRestrictTargetAddressing(fineLevelIndex);

    const labelgpuListList& patchFineToCoarseTargetStart =
        agglomeration_.patchFaceRestrictTargetStartAddressing(fineLevelIndex);

    const labelList& nPatchFaces =
        agglomeration_.nPatchFaces(fineLevelIndex);


    // Add the coarse level
    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            const GAMGInterface& coarseInterface =
                refCast<const GAMGInterface>
                (
                    coarseMeshInterfaces[inti]
                );

            coarsePrimInterfaces.set
            (
                inti,
                GAMGInterfaceField::New
                (
                    coarseInterface,
                    fineInterfaces[inti]
                ).ptr()
            );
            coarseInterfaces.set
            (
                inti,
                &coarsePrimInterfaces[inti]
            );

            const labelgpuList& faceRestrictSortAddressing = patchFineToCoarseSort[inti];
            const labelgpuList& faceRestrictTargetAddressing = patchFineToCoarseTarget[inti];
            const labelgpuList& faceRestrictTargetStartAddressing = patchFineToCoarseTargetStart[inti];

            coarseInterfaceBouCoeffs.set
            (
                inti,
                new scalargpuField(nPatchFaces[inti], 0.0)
            );
            agglomeration_.restrictField
            (
                coarseInterfaceBouCoeffs[inti],
                fineInterfaceBouCoeffs[inti],
                faceRestrictSortAddressing,
                faceRestrictTargetAddressing,
                faceRestrictTargetStartAddressing
            );

            coarseInterfaceIntCoeffs.set
            (
                inti,
                new scalargpuField(nPatchFaces[inti], 0.0)
            );
            agglomeration_.restrictField
            (
                coarseInterfaceIntCoeffs[inti],
                fineInterfaceIntCoeffs[inti],
                faceRestrictSortAddressing,
                faceRestrictTargetAddressing,
                faceRestrictTargetStartAddressing
            );
        }
    }
}


// Gather matrices.
// Note: matrices get constructed with dummy mesh
void Foam::GAMGSolver::gatherMatrices
(
    const labelList& procIDs,
    const lduMesh& dummyMesh,
    const label meshComm,

    const lduMatrix& mat,
    const FieldField<gpuField, scalar>& interfaceBouCoeffs,
    const FieldField<gpuField, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,

    PtrList<lduMatrix>& otherMats,
    PtrList<FieldField<gpuField, scalar> >& otherBouCoeffs,
    PtrList<FieldField<gpuField, scalar> >& otherIntCoeffs,
    List<boolList>& otherTransforms,
    List<List<int> >& otherRanks
) const
{
    if (debug)
    {
        Pout<< "GAMGSolver::gatherMatrices :"
            << " collecting matrices from procs:" << procIDs
            << " using comm:" << meshComm << endl;
    }

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Master.
        otherMats.setSize(procIDs.size()-1);
        otherBouCoeffs.setSize(procIDs.size()-1);
        otherIntCoeffs.setSize(procIDs.size()-1);
        otherTransforms.setSize(procIDs.size()-1);
        otherRanks.setSize(procIDs.size()-1);

        for (label procI = 1; procI < procIDs.size(); procI++)
        {
            label otherI = procI-1;

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[procI],
                0,          // bufSize
                Pstream::msgType(),
                meshComm
            );

            otherMats.set(otherI, new lduMatrix(dummyMesh, fromSlave));

            // Receive number of/valid interfaces
            boolList& procTransforms = otherTransforms[otherI];
            List<int>& procRanks = otherRanks[otherI];

            fromSlave >> procTransforms;
            fromSlave >> procRanks;

            // Size coefficients
            otherBouCoeffs.set
            (
                otherI,
                new FieldField<gpuField, scalar>(procRanks.size())
            );
            otherIntCoeffs.set
            (
                otherI,
                new FieldField<gpuField, scalar>(procRanks.size())
            );
            forAll(procRanks, intI)
            {
                if (procRanks[intI] != -1)
                {
                    otherBouCoeffs[otherI].set
                    (
                        intI,
                        new scalargpuField(fromSlave)
                    );
                    otherIntCoeffs[otherI].set
                    (
                        intI,
                        new scalargpuField(fromSlave)
                    );
                }
            }
        }
    }
    else
    {
        // Send to master

        // Count valid interfaces
        boolList procTransforms(interfaceBouCoeffs.size(), false);
        List<int> procRanks(interfaceBouCoeffs.size(), -1);
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterfaceField& interface =
                    refCast<const processorLduInterfaceField>
                    (
                        interfaces[intI]
                    );

                procTransforms[intI] = interface.doTransform();
                procRanks[intI] = interface.rank();
            }
        }

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            meshComm
        );

        toMaster << mat << procTransforms << procRanks;
        forAll(procRanks, intI)
        {
            if (procRanks[intI] != -1)
            {
                toMaster
                    << interfaceBouCoeffs[intI]
                    << interfaceIntCoeffs[intI];
            }
        }
    }
}

namespace Foam
{

struct luProcAgglomerateFunctor
{
    scalar* lu;
    luProcAgglomerateFunctor(scalar* lu_):lu(lu_){}

    __HOST____DEVICE__
    void operator()(const thrust::tuple<scalar,scalar,scalar>& t)
    {
        label map = thrust::get<0>(t);
        if(map>=0)
        {
            lu[map] = -thrust::get<1>(t);
        }
        else
        {
            lu[-map-1] = -thrust::get<2>(t);
        }
    }
};

}

void Foam::GAMGSolver::procAgglomerateMatrix
(
    // Agglomeration information
    const labelList& procAgglomMap,
    const List<int>& agglomProcIDs,

    const label levelI,

    // Resulting matrix
    autoPtr<lduMatrix>& allMatrixPtr,
    FieldField<gpuField, scalar>& allInterfaceBouCoeffs,
    FieldField<gpuField, scalar>& allInterfaceIntCoeffs,
    PtrList<lduInterfaceField>& allPrimitiveInterfaces,
    lduInterfaceFieldPtrsList& allInterfaces
) const
{
    const lduMatrix& coarsestMatrix = matrixLevels_[levelI];
    const lduInterfaceFieldPtrsList& coarsestInterfaces =
        interfaceLevels_[levelI];
    const FieldField<gpuField, scalar>& coarsestBouCoeffs =
        interfaceLevelsBouCoeffs_[levelI];
    const FieldField<gpuField, scalar>& coarsestIntCoeffs =
        interfaceLevelsIntCoeffs_[levelI];
    const lduMesh& coarsestMesh = coarsestMatrix.mesh();


    label coarseComm = coarsestMesh.comm();

    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = coarseComm;



    // Gather all matrix coefficients onto agglomProcIDs[0]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PtrList<lduMatrix> otherMats;
    PtrList<FieldField<gpuField, scalar> > otherBouCoeffs;
    PtrList<FieldField<gpuField, scalar> > otherIntCoeffs;
    List<boolList> otherTransforms;
    List<List<int> > otherRanks;
    gatherMatrices
    (
        agglomProcIDs,
        coarsestMesh,
        coarseComm,

        coarsestMatrix,
        coarsestBouCoeffs,
        coarsestIntCoeffs,
        coarsestInterfaces,

        otherMats,
        otherBouCoeffs,
        otherIntCoeffs,
        otherTransforms,
        otherRanks
    );


    if (Pstream::myProcNo(coarseComm) == agglomProcIDs[0])
    {
        // Agglomerate all matrix
        // ~~~~~~~~~~~~~~~~~~~~~~

        //Pout<< "Own matrix:" << coarsestMatrix.info() << endl;
        //
        //forAll(otherMats, i)
        //{
        //    Pout<< "** otherMats " << i << " "
        //        << otherMats[i].info()
        //        << endl;
        //}
        //Pout<< endl;


        const lduMesh& allMesh = agglomeration_.meshLevel(levelI+1);
        const labelList& cellOffsets = agglomeration_.cellOffsets(levelI+1);
        const labelgpuListList& faceMap = agglomeration_.faceMap(levelI+1);
        const labelListList& boundaryMap = agglomeration_.boundaryMap(levelI+1);
        const labelgpuListListList& boundaryFaceMap =
            agglomeration_.boundaryFaceMap(levelI+1);

        allMatrixPtr.reset(new lduMatrix(allMesh));
        lduMatrix& allMatrix = allMatrixPtr();

        if (coarsestMatrix.hasDiag())
        {
            scalargpuField& allDiag = allMatrix.diag();
            gpuList<scalar>
            (
                const_cast<const scalargpuField&>(allDiag),
                coarsestMatrix.diag().size()
            )
            = coarsestMatrix.diag();
            forAll(otherMats, i)
            {
                gpuList<scalar>
                (
                    allDiag,
                    otherMats[i].diag().size(),
                    cellOffsets[i+1]
                )
                = otherMats[i].diag();
            }
        }
        if (coarsestMatrix.hasLower())
        {
            scalargpuField& allLower = allMatrix.lower();
            gpuIndirectList<scalar>
            (
                allLower,
                faceMap[0]
            ) = coarsestMatrix.lower();
            forAll(otherMats, i)
            {
                gpuIndirectList<scalar>
                (
                    allLower,
                    faceMap[i+1]
                ) = otherMats[i].lower();
            }
        }
        if (coarsestMatrix.hasUpper())
        {
            scalargpuField& allUpper = allMatrix.upper();
            gpuIndirectList<scalar>
            (
                allUpper,
                faceMap[0]
            ) = coarsestMatrix.upper();
            forAll(otherMats, i)
            {
                gpuIndirectList<scalar>
                (
                    allUpper,
                    faceMap[i+1]
                ) = otherMats[i].upper();
            }
        }


        // Agglomerate interface fields and coefficients
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        lduInterfacePtrsList allMeshInterfaces = allMesh.interfaces();

        allInterfaceBouCoeffs.setSize(allMeshInterfaces.size());
        allInterfaceIntCoeffs.setSize(allMeshInterfaces.size());
        allPrimitiveInterfaces.setSize(allMeshInterfaces.size());
        allInterfaces.setSize(allMeshInterfaces.size());

        forAll(allMeshInterfaces, intI)
        {
            const lduInterface& patch = allMeshInterfaces[intI];
            label size = patch.faceCells().size();

            allInterfaceBouCoeffs.set(intI, new scalargpuField(size));
            allInterfaceIntCoeffs.set(intI, new scalargpuField(size));
        }

        labelList nBounFaces(allMeshInterfaces.size());
        forAll(boundaryMap, procI)
        {
            const FieldField<gpuField, scalar>& procBouCoeffs
            (
                (procI == 0)
              ? coarsestBouCoeffs
              : otherBouCoeffs[procI-1]
            );
            const FieldField<gpuField, scalar>& procIntCoeffs
            (
                (procI == 0)
              ? coarsestIntCoeffs
              : otherIntCoeffs[procI-1]
            );

            const labelList& bMap = boundaryMap[procI];
            forAll(bMap, procIntI)
            {
                label allIntI = bMap[procIntI];

                if (allIntI != -1)
                {
                    // So this boundary has been preserved. Copy
                    // data across.

                    if (!allInterfaces.set(allIntI))
                    {
                        // Construct lduInterfaceField

                        bool doTransform = false;
                        int rank = -1;
                        if (procI == 0)
                        {
                            const processorGAMGInterfaceField& procInt =
                                refCast
                                <
                                    const processorGAMGInterfaceField
                                >
                                (
                                    coarsestInterfaces[procIntI]
                                );
                            doTransform = procInt.doTransform();
                            rank = procInt.rank();
                        }
                        else
                        {
                            doTransform =
                                otherTransforms[procI-1][procIntI];
                            rank = otherRanks[procI-1][procIntI];
                        }

                        allPrimitiveInterfaces.set
                        (
                            allIntI,
                            GAMGInterfaceField::New
                            (
                                refCast<const GAMGInterface>
                                (
                                    allMeshInterfaces[allIntI]
                                ),
                                doTransform,
                                rank
                            ).ptr()
                        );
                        allInterfaces.set
                        (
                            allIntI,
                            &allPrimitiveInterfaces[allIntI]
                        );
                    }


                    // Map data from processor to complete mesh

                    scalargpuField& allBou = allInterfaceBouCoeffs[allIntI];
                    scalargpuField& allInt = allInterfaceIntCoeffs[allIntI];

                    const labelgpuList& map = boundaryFaceMap[procI][procIntI];

                    const scalargpuField& procBou = procBouCoeffs[procIntI];
                    const scalargpuField& procInt = procIntCoeffs[procIntI];

                    thrust::copy
                    (
                        procBou.begin(),
                        procBou.end(),
                        thrust::make_permutation_iterator
                        (
                            allBou.begin(),
                            map.begin()
                        )
                    );

                    thrust::copy
                    (
                        procInt.begin(),
                        procInt.end(),
                        thrust::make_permutation_iterator(
                            allInt.begin(),
                            map.begin()
                        )
                    );
                }
                else if (procBouCoeffs.set(procIntI))
                {
                    // Boundary has become internal face

                    const labelgpuList& map = boundaryFaceMap[procI][procIntI];
                    const scalargpuField& procBou = procBouCoeffs[procIntI];
                    const scalargpuField& procInt = procIntCoeffs[procIntI];

                    if(coarsestMatrix.hasUpper())
                    {
                        thrust::for_each
                        (
                            thrust::make_zip_iterator(thrust::make_tuple
                            (
                                map.begin(),
                                procBou.begin(),
                                procInt.begin()
                            )),
                            thrust::make_zip_iterator(thrust::make_tuple
                            (
                                map.end(),
                                procBou.end(),
                                procInt.end()
                            )),
                            luProcAgglomerateFunctor(allMatrix.upper().data())
                        );

                    }

                    if(coarsestMatrix.hasLower())
                    {
                        thrust::for_each
                        (
                            thrust::make_zip_iterator(thrust::make_tuple
                            (
                                map.begin(),
                                procInt.begin(),
                                procBou.begin()
                            )),
                            thrust::make_zip_iterator(thrust::make_tuple
                            (
                                map.end(),
                                procInt.end(),
                                procBou.end()
                            )),
                            luProcAgglomerateFunctor(allMatrix.lower().data())
                        );

                    }
                }
            }
        }

        //Pout<< "** Assembled allMatrix:" << allMatrix.info() << endl;
        //
        //forAll(allInterfaces, intI)
        //{
        //    if (allInterfaces.set(intI))
        //    {
        //        Pout<< "    patch:" << intI
        //            << " type:" << allInterfaces[intI].type()
        //            << " size:"
        //            << allInterfaces[intI].interface().
        //                faceCells().size()
        //            << endl;
        //
        //        //const scalarField& bouCoeffs = allInterfaceBouCoeffs[intI];
        //        //const scalarField& intCoeffs = allInterfaceIntCoeffs[intI];
        //        //forAll(bouCoeffs, faceI)
        //        //{
        //        //    Pout<< "        " << faceI
        //        //        << "\tbou:" << bouCoeffs[faceI]
        //        //        << "\tint:" << intCoeffs[faceI]
        //        //        << endl;
        //        //}
        //    }
        //}
    }
    UPstream::warnComm = oldWarn;
}


void Foam::GAMGSolver::procAgglomerateMatrix
(
    const labelList& procAgglomMap,
    const List<int>& agglomProcIDs,

    const label levelI
)
{
    autoPtr<lduMatrix> allMatrixPtr;
    autoPtr<FieldField<gpuField, scalar> > allInterfaceBouCoeffs
    (
        new FieldField<gpuField, scalar>(0)
    );
    autoPtr<FieldField<gpuField, scalar> > allInterfaceIntCoeffs
    (
        new FieldField<gpuField, scalar>(0)
    );
    autoPtr<PtrList<lduInterfaceField> > allPrimitiveInterfaces
    (
        new PtrList<lduInterfaceField>(0)
    );
    autoPtr<lduInterfaceFieldPtrsList> allInterfaces
    (
        new lduInterfaceFieldPtrsList(0)
    );

    procAgglomerateMatrix
    (
        // Agglomeration information
        procAgglomMap,
        agglomProcIDs,

        levelI,

        // Resulting matrix
        allMatrixPtr,
        allInterfaceBouCoeffs(),
        allInterfaceIntCoeffs(),
        allPrimitiveInterfaces(),
        allInterfaces()
    );

    matrixLevels_.set(levelI, allMatrixPtr);
    interfaceLevelsBouCoeffs_.set(levelI, allInterfaceBouCoeffs);
    interfaceLevelsIntCoeffs_.set(levelI, allInterfaceIntCoeffs);
    primitiveInterfaceLevels_.set(levelI, allPrimitiveInterfaces);
    interfaceLevels_.set(levelI, allInterfaces);
}


// ************************************************************************* //
