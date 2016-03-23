/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::radiation::wsggmAbsorptionEmission

Description
    wsggmAbsorptionEmission radiation absorption and emission coefficients
    for continuous phase

    Dmitry Lysenko (c) 2014 	
    Rev. 17 August 2014 	    

SourceFiles
    wsggmAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wsggmAbsorptionEmission_H
#define wsggmAbsorptionEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "fluidThermo.H"
#include "fvCFD.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
   Rev. 21 August 2014
   Weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission

   Rev. 6 September 2014
   Added a function to calculate the mean beam lengths for optically thin media,
   L0 = 4 V/A (where A is the entire area bounding the volume V), which are often called geometric mean beam lengths,
   based on the work by Dunkle [1964].

   Dunkle, R. V.: "Geometric mean beam lengths for radiant heat transfer calculations,"
   ASME Journal of Heat Transfer, vol. 86, no. 1, pp. 75-80, 1964

   Rev. 27 September 2014
   Add a possibility to swith between manual (setup in config dict) and auto calculation of the mean beam path.
   Actual for 2D planar cases;
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class greyMeanAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggmAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- absorption coefficitent for individual grey gases
        scalarList emissivityCoeffs_;

        //- fitting Factors
        scalarListList fittingFactors_;

        //characteristic path length / [m] from dict
        dimensionedScalar pathLength_;

        //- SLG thermo package
        const fluidThermo& thermo_;
        
        dimensionedScalar sector_;  // [deg] from dict


        // atomic weights
        const scalar _wCO2;
        const scalar _wH2O;

        // calc mean beam path length:
        // true: automatically
        // false: manual / read the value from the coeffsDdict
        bool meanBeamPathAutoCalcMode_;  // from dict


public:

    //- Runtime type information
    TypeName("wsggmAbsorptionEmission");


    // Constructors

        //- Construct from components
        wsggmAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggmAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
