
Info<< "Creating combustion model\n" << endl;


autoPtr<combustionModels::rhoCombustionModel> combustion
(
    combustionModels::rhoCombustionModel::New(mesh)
);

rhoReactionThermo& thermo = combustion->thermo(); /*  hsReactionThermo& thermo = chemistry.thermo(); */ 
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition(); /*  --- same as OF 2.1.1 --- PtrList<volScalarField>& Y = composition.Y(); */ 
PtrList<volScalarField>& Y = composition.Y();  /*  --- same as OF 2.1.1 --- PtrList<volScalarField>& Y = composition.Y(); */ 


/*************************************************************************************************************/

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
	IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField Cp
(
    IOobject
    (
        "thermo:Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);

volScalarField Cv
(
    IOobject
    (
        "thermo:Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cv()
);
volScalarField kappa
(
    IOobject
    (
        "thermo:kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);

volScalarField mu
(
    IOobject
    (
        "thermo:mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
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


volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
volScalarField& he = thermo.he();	/* instead of  volScalarField& hs = thermo.hs();  */
const volScalarField& T = thermo.T();


#include "compressibleCreatePhi.H"


Info << "Creating turbulence model\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

combustion->setTurbulence(turbulence());

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he()); /* instead of fields.add(hs); */ 

Info<< "Creating ignition source\n" << endl;
IOdictionary ignitionProperties
(
       IOobject
       (
           "ignitionProperties",
           runTime.constant(),
           mesh,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       )
);

ignition ign(ignitionProperties, runTime, mesh);


Info<< "Creating additional terms for EDC-LES\n" << endl;


volScalarField tauK
(
    IOobject
    (
        "tauK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tauK", dimensionSet(0,0,1,0,0,0,0), 0.0)
);


volScalarField tauDelta
(
    IOobject
    (
        "tauDelta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tauDelta", dimensionSet(0,0,1,0,0,0,0), 0.0)
);

volScalarField tauStar
(
    IOobject
    (
        "tauStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tauStar", dimensionSet(0,0,1,0,0,0,0), 0.0)
);


const dimensionedScalar tauMin_( "tauStarMin_", dimensionSet(0,0,1,0,0,0,0), SMALL ); 
const dimensionedScalar eSgsMin_( "eSgsMin_",   dimensionSet(0,2,-3,0,0,0,0), SMALL ); 

volScalarField cDelta
(
    IOobject
    (
        "cDelta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cDelta", dimensionSet(0,1,0,0,0,0,0), 0.0)
);


cDelta.internalField() = pow ( mesh.V(), 1.0/3.0) ;


volScalarField eSgs
(
    IOobject
    (
        "eSgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eSgs", dimensionSet(0,2,-3,0,0,0,0), 1e-3)
);

volScalarField kSgs
(
    IOobject
    (
        "kSgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kSgs", dimensionSet(0,2,-2,0,0,0,0), 1e-3)
);


volScalarField gammaL
(
    IOobject
    (
        "gammaL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gammaL", dimensionSet(0,0,0,0,0,0,0), 0.0)
);


volScalarField ReTauLES_
(
	    IOobject
	    (
	        "ReTauLES_",
	        runTime.timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("ReTauLES_", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

