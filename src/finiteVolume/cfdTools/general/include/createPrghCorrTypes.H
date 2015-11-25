    Foam::wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        Foam::zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p_rgh.boundaryField().size(); i++)
    {
        if (p_rgh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = Foam::fixedValueFvPatchScalarField::typeName;
        }
    }
