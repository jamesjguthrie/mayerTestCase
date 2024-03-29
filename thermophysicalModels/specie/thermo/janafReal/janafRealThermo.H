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
    Foam::janafRealThermo

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    janafRealThermoI.H
    janafRealThermo.C

\*---------------------------------------------------------------------------*/

#ifndef janafRealThermo_H
#define janafRealThermo_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class janafRealThermo;

template<class EquationOfState>
inline janafRealThermo<EquationOfState> operator+
(
    const janafRealThermo<EquationOfState>&,
    const janafRealThermo<EquationOfState>&
);

template<class EquationOfState>
inline janafRealThermo<EquationOfState> operator-
(
    const janafRealThermo<EquationOfState>&,
    const janafRealThermo<EquationOfState>&
);

template<class EquationOfState>
inline janafRealThermo<EquationOfState> operator*
(
    const scalar,
    const janafRealThermo<EquationOfState>&
);

template<class EquationOfState>
inline janafRealThermo<EquationOfState> operator==
(
    const janafRealThermo<EquationOfState>&,
    const janafRealThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const janafRealThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janafRealThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class janafRealThermo
:
    public EquationOfState
{

public:

    // Public data

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline janafRealThermo
        (
            const EquationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highCpCoeffs,
            const coeffArray& lowCpCoeffs
        );

        //- Construct from Istream
        janafRealThermo(Istream&);

        //- Construct from dictionary
        janafRealThermo(const dictionary& dict);

        //- Construct as a named copy
        inline janafRealThermo(const word&, const janafRealThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "janafReal<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the common temperature
            inline scalar Tcommon() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowCpCoeffs() const;


        // Fundamental properties

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


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const janafRealThermo&);
        inline void operator-=(const janafRealThermo&);


    // Friend operators

        friend janafRealThermo operator+ <EquationOfState>
        (
            const janafRealThermo&,
            const janafRealThermo&
        );

        friend janafRealThermo operator- <EquationOfState>
        (
            const janafRealThermo&,
            const janafRealThermo&
        );

        friend janafRealThermo operator* <EquationOfState>
        (
            const scalar,
            const janafRealThermo&
        );

        friend janafRealThermo operator== <EquationOfState>
        (
            const janafRealThermo&,
            const janafRealThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const janafRealThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janafRealThermoI.H"

#ifdef NoRepository
#   include "janafRealThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
