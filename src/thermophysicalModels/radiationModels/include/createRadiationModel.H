    autoPtr<Foam::radiation::radiationModel> radiation
    (
        Foam::radiation::radiationModel::New(thermo.T())
    );
