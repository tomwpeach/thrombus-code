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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar dTh
    (
        transportProperties.lookup("dTh")
    );

    dimensionedScalar ThMin
    (
        transportProperties.lookup("ThMin")
    );
    
    dimensionedScalar ThCrit
    (
        transportProperties.lookup("ThCrit")
    );

    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );

    dimensionedScalar e
    (
        transportProperties.lookup("e")
    );
    
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Th\n" << endl;
    volScalarField Th
    (
        IOobject
        (
            "Th",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    

    Info<< "Reading field s\n" << endl;
    volScalarField s
    (
        IOobject
        (
            "s",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading field m\n" << endl;
    volScalarField m
    (
        IOobject
        (
            "m",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field mu\n" << endl;
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field wallGradU\n" << endl;
    volVectorField wallGradU
            (
                IOobject
                (
                    "wallGradU",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "wallGradU",
                    U.dimensions()/dimLength,
                    vector::zero
                )
	     );


    Info<< "Reading field StrainRateMag\n" << endl;
    volScalarField StrainRateMag
    (
        IOobject
        (
            "StrainRateMag",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );		
    
#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
