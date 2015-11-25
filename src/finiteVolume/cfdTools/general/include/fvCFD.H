#ifndef fvCFD_H
#define fvCFD_H

#include "parRun.H"

#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "constants.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

#ifndef namespaceFoam
#define namespaceFoam

    #include "accessFoamNamespace.H"

    using Foam::fvMesh;
    using Foam::volMesh;
    using Foam::fvMatrix;
    using Foam::fvScalarMatrix;
    using Foam::fvVectorMatrix;
    using Foam::solve;

    using Foam::fvSolution;

    using Foam::volScalarField;
    using Foam::volVectorField;
    using Foam::volTensorField;
    using Foam::surfaceScalarField;
    using Foam::surfaceVectorField;
    using Foam::fvPatchField;
    using Foam::fvsPatchScalarField;
    using Foam::volSymmTensorField;
    using Foam::uniformDimensionedVectorField;

    namespace fvm = Foam::fvm;
    namespace fvc = Foam::fvc;

#endif

#endif
