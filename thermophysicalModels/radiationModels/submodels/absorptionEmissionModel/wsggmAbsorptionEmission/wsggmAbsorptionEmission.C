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

#include "wsggmAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::wsggmAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    emissivityCoeffs_(coeffsDict_.lookup("emissivityCoeffs")),
    fittingFactors_(coeffsDict_.lookup("fittingFactors")),
    pathLength_(coeffsDict_.lookup("pathLength")),
    thermo_(mesh.lookupObject<fluidThermo>("thermophysicalProperties")),
    _wCO2(44),
    _wH2O(18),
    sector_(coeffsDict_.lookup("sector"))
   
{
    if (!isA<basicMultiComponentMixture>(thermo_))
    {
        FatalErrorIn
        (
            "radiation::wsggmAbsorptionEmission::wsggmAbsorptionEmission"
            "("
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    const scalar s_ = sector_.value();


    if (s_ <= 0.0 || s_>360.0)
       {
           FatalErrorIn
           (
               "radiation::wsggmAbsorptionEmission::wsggmAbsorptionEmission"
               "("
                   "const dictionary&, "
                   "const fvMesh&"
               ")"
           )   << "Sector should be specified in a range:   0.0 < s << 360.0 "  << abort(FatalError);
       }


    // calc mean beam path length:
     // true: automatically
     // false: manual / read the value from the coeffsDdict

    meanBeamPathAutoCalcMode_       = readBool(coeffsDict_.lookup("meanBeamPathAutoCalcMode"));


    // calculate mean beam path length automatically

    if   (meanBeamPathAutoCalcMode_)
    {

       scalar dA = 0.0;
       scalar dAall = 0.0;
       scalar dAw = 0.0;
       scalar dAs = 0.0;

       scalar dV = 0.0;

       const fvPatchList& boundaries = mesh.boundary();


       dV  = sum(mesh.V().field());
       reduce(dV, sumOp<scalar>());







       for (int patchI=0; patchI < boundaries.size(); patchI++)
       {

           if ( boundaries[patchI].type() != "processor" )
                 dAall += sum(mesh.magSf().boundaryField()[patchI]);

           if ( boundaries[patchI].type() == "wedge" )
                 dAw += sum(mesh.magSf().boundaryField()[patchI]);

           if ( boundaries[patchI].type() == "symmetryPlane" )
                 dAs += sum(mesh.magSf().boundaryField()[patchI]);

       }
       reduce(dAall, sumOp<scalar>());
       reduce(dAw, sumOp<scalar>());
       reduce(dAs, sumOp<scalar>());


       dA = dAall - dAw - dAs;



       // recalculate dV and dA from sectror to full 3d values
       if (s_ != 360.0)
       {

             dV = 360.0/s_ * dV;
             dA = 360.0/s_ * dA;
       }


       dimensionedScalar domainArea("dA",dimensionSet(0,2,0,0,0,0,0),dA);
       dimensionedScalar domainVolume("dV",dimensionSet(0,3,0,0,0,0,0),dV);
       dimensionedScalar factor("4",dimensionSet(0,0,0,0,0,0,0),scalar(4.0));


       pathLength_ = factor * domainVolume / domainArea;
       Info << "Domain volume (m3)  :  " << domainVolume.value()  <<  endl;
       Info << "Domain area (m2)    :  " << domainArea.value()  <<  endl;

     } // end if



     Info << "Mean beam length (m) for WSGGM was set to " << pathLength_.value()  <<  endl;










    /*
    const fvPatchList& boundaries = mesh.boundary();

    scalar dA = 0.0;
    forAll(boundaries, patchI)
    { 
       // we are not calculating area of the wedges   
       if (boundaries[patchI].type() != "wedge")		
	    dA += gSum(mesh.magSf().boundaryField()[patchI]);
    } // end for boundaries 


    scalar dV = 0.0;
    forAll(mesh.V(), cellI)
    { 	
	dV += mesh.V()[cellI];
    }   



     
    // recalculate dV and dA from sectror to full 3d values
    if (s_ != 360.0)
    {
	
    	dV = 360.0/s_ * dV;
	dA = 360.0/s_ * dA;	     
    }

	
    dimensionedScalar domainArea("dA",dimensionSet(0,2,0,0,0,0,0),dA); 
    dimensionedScalar domainVolume("dV",dimensionSet(0,3,0,0,0,0,0),dV);
    dimensionedScalar factor("4",dimensionSet(0,0,0,0,0,0,0),scalar(4.0));   

    
    pathLength_ = factor * domainVolume / domainArea; 
	
    if (s_ != 360) 
	Info << "Domain type is a sector of " << sector_.value() <<  " deg" <<  endl;	     
    else
	Info << "Domain type is a full 360 deg domain" <<  endl;	     	

    Info << "Domain volume (m3)  :  " << domainVolume.value()  <<  endl;	 
    Info << "Domain area (m2)    :  " << domainArea.value()  <<  endl;	 
    Info << "Mean beam length (m) for WSGGM was set to " << pathLength_.value()  <<  endl;
    */
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::~wsggmAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::aCont(const label bandI) const
{
    const basicMultiComponentMixture& mixture =
        dynamic_cast<const basicMultiComponentMixture&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    const label indexCO2= mixture.species()["CO2"];
    const label indexH2O= mixture.species()["H2O"];

    //const scalar wCO2=44;
    //const scalar wH2O=18;


   
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );


    volScalarField emissivity(        
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    volScalarField pressurePathLength(        
            IOobject
            (
                "pressurePathLength",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimPressure*dimLength,0.0)
        );
    volScalarField weightingFactor(        
            IOobject
            (
                "weightingFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );

    volScalarField wMean(        
            IOobject
            (
                "wMean",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0),0.0) // kg/kmol
        );


    forAll (mixture.Y(), s)
    {
        wMean += mixture.Y(s)/mixture.W(s);
    }

  

    wMean=1/wMean;


   

    pressurePathLength = wMean*p*(mixture.Y(indexCO2)/mixture.W(indexCO2) + mixture.Y(indexH2O)/mixture.W(indexH2O))*pathLength_;

    volScalarField limitedDimlessTemperature = min( T /dimensionedScalar("unity",dimTemperature, 1.0), 2400.0);

    forAll(emissivityCoeffs_,i)
    {
        weightingFactor = 0.0;
        for(label j=1; j<fittingFactors_.size(); j++)
        {
            weightingFactor += fittingFactors_[i][j]*pow(limitedDimlessTemperature, (j-1));
        }
        emissivity += weightingFactor * (1 - exp( (-1) *emissivityCoeffs_[i] * pathLength_.value()));
    }

    emissivity = min(emissivity,0.9999);

    ta()= (-1)* log(1-emissivity) / pathLength_;



    ta().correctBoundaryConditions();
    return ta;


  
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::eCont(const label bandI) const
{
   return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::ECont(const label bandI) const
{

  
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
   
    return E;
 
}


// ************************************************************************* //
