    Foam::wordList pcorrTypes
    (
        p.boundaryField().size(),
        Foam::zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = Foam::fixedValueFvPatchScalarField::typeName;
        }
    }
