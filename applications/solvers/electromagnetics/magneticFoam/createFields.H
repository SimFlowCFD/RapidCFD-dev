    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    List<magnet> magnets(transportProperties.lookup("magnets"));

    surfaceScalarField murf
    (
        IOobject
        (
            "murf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        1
    );

    surfaceScalarField Mrf
    (
        IOobject
        (
            "Mrf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("Mr", dimensionSet(0, 1, 0, 0, 0, 1, 0), 0)
    );

    forAll(magnets, i)
    {
        label magnetZonei = mesh.faceZones().findZoneID(magnets[i].name());

        if (magnetZonei == -1)
        {
            FatalIOErrorIn(args.executable().c_str(), transportProperties)
                << "Cannot find faceZone for magnet " << magnets[i].name()
                << exit(FatalIOError);
        }

        const labelgpuList& faces = mesh.faceZones()[magnetZonei].getList();

        const scalar muri = magnets[i].mur();
        const scalar Mri = magnets[i].Mr().value();
        const vector& orientationi = magnets[i].orientation();

        const surfaceVectorField& Sf = mesh.Sf();

        thrust::fill
        (
            thrust::make_permutation_iterator
            (
                murf.getField().begin(),
                faces.begin()
            ),
            thrust::make_permutation_iterator
            (
                murf.getField().begin(),
                faces.end()
            ),
            muri
        );

        thrust::transform
        (
            thrust::make_constant_iterator(Mri),
            thrust::make_constant_iterator(Mri)+faces.size(),
            thrust::make_transform_iterator
            (
                thrust::make_permutation_iterator
                (
                    Sf.getField().begin(),
                    faces.begin()
                ),
                Foam::dotOperatorSFFunctor<vector,vector,scalar>(orientationi)
            ),
            thrust::make_permutation_iterator
            (
                Mrf.getField().begin(),
                faces.begin()
            ),
            Foam::multiplyOperatorFunctor<scalar,scalar,scalar>()
        );
    }
