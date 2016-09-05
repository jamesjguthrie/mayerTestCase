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
    Foam::CpFromEoSThermo

SourceFiles
    CpFromEoSThermoI.H
    CpFromEoSThermo.C

\*---------------------------------------------------------------------------*/

/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Calculate cp, cv, etc. directly    
*/

#ifndef CpFromEoSThermo_H
#define CpFromEoSThermo_H

#include "scalar.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState, int PolySize>
class CpFromEoSThermo;

template<class EquationOfState, int PolySize>
inline CpFromEoSThermo<EquationOfState, PolySize> operator+
(
    const CpFromEoSThermo<EquationOfState, PolySize>&,
    const CpFromEoSThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline CpFromEoSThermo<EquationOfState, PolySize> operator-
(
    const CpFromEoSThermo<EquationOfState, PolySize>&,
    const CpFromEoSThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline CpFromEoSThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const CpFromEoSThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline CpFromEoSThermo<EquationOfState, PolySize> operator==
(
    const CpFromEoSThermo<EquationOfState, PolySize>&,
    const CpFromEoSThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const CpFromEoSThermo<EquationOfState, PolySize>&
);


/*---------------------------------------------------------------------------*\
                      Class CpFromEoSThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize=8>
class CpFromEoSThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;

        //- Standard entropy
        //  Note: input in [J/kg/K], but internally uses [J/kmol/K]
        scalar Sf_;

        //- Specific heat at constant pressure polynomial coeffs [J/(kg.K)]
        Polynomial<PolySize> CpCoeffs_;

	mutable scalar counter;

	//- Enthalpy polynomial coeffs - derived from cp [J/kg]
        //  NOTE: relative to Tstd
        typename Polynomial<PolySize>::intPolyType hCoeffs_;

        //- Entropy - derived from Cp [J/(kg.K)] - relative to Tstd
        Polynomial<PolySize> sCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline CpFromEoSThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
	    const Polynomial<PolySize>& CpCoeffs,
            const typename Polynomial<PolySize>::intPolyType& hCoeffs,
            const Polynomial<PolySize>& sCoeffs
	);


public:

    // Constructors

        //- Construct from Istream
        CpFromEoSThermo(Istream& is);

        //- Construct from dictionary
        CpFromEoSThermo(const dictionary& dict);

        //- Construct as copy
        inline CpFromEoSThermo(const CpFromEoSThermo&);

        //- Construct as a named copy
        inline CpFromEoSThermo(const word&, const CpFromEoSThermo&);

    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CpFromEoS<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties
        
	    //- perfect Gas internal Energy  [J/kmol]
            inline scalar e0(const scalar T) const;
        
	    //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

	    inline scalar cv(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;
	    
            //- Internal Energy [J/kmol]
            inline scalar e(const scalar rho, const scalar T) const;

        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline CpFromEoSThermo& operator=(const CpFromEoSThermo&);
        inline void operator+=(const CpFromEoSThermo&);
        inline void operator-=(const CpFromEoSThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend CpFromEoSThermo operator+ <EquationOfState, PolySize>
        (
            const CpFromEoSThermo&,
            const CpFromEoSThermo&
        );

        friend CpFromEoSThermo operator- <EquationOfState, PolySize>
        (
            const CpFromEoSThermo&,
            const CpFromEoSThermo&
        );

        friend CpFromEoSThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const CpFromEoSThermo&
        );

        friend CpFromEoSThermo operator== <EquationOfState, PolySize>
        (
            const CpFromEoSThermo&,
            const CpFromEoSThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const CpFromEoSThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CpFromEoSThermoI.H"

#ifdef NoRepository
#   include "CpFromEoSThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
