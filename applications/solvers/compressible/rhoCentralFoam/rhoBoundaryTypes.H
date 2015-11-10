const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
Foam::wordList rhoBoundaryTypes = pbf.types();

forAll(rhoBoundaryTypes, patchi)
{
    if (rhoBoundaryTypes[patchi] == "waveTransmissive")
    {
        rhoBoundaryTypes[patchi] = Foam::zeroGradientFvPatchScalarField::typeName;
    }
    else if (pbf[patchi].fixesValue())
    {
        rhoBoundaryTypes[patchi] = Foam::fixedRhoFvPatchScalarField::typeName;
    }
}
