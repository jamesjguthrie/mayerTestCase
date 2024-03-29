Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::rhoCombustionModel> reaction
(
    combustionModels::rhoCombustionModel::New(mesh)
);

rhoReactionThermo& thermo = reaction->thermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

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

Cp.write();

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

Cv.write();

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
const volScalarField& T = thermo.T();

/*volScalarField psiByDp
(
    IOobject
    (
        "psiByDp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (psi * p)
);*/


#include "compressibleCreatePhi.H"


Info << "Creating turbulence model.\n" << nl;
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

// Set the turbulence into the reaction model
reaction->setTurbulence(turbulence());


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);
