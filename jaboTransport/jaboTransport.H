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
    Foam::jaboTransport

Description
    Transport package using jabo functions for mu and kappa

SourceFiles
    jaboTransportI.H
    jaboTransport.C

\*---------------------------------------------------------------------------*/

#ifndef jaboTransport_H
#define jaboTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class jaboTransport;

template<class Thermo, int PolySize>
inline jaboTransport<Thermo, PolySize> operator+
(
    const jaboTransport<Thermo, PolySize>&,
    const jaboTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline jaboTransport<Thermo, PolySize> operator-
(
    const jaboTransport<Thermo, PolySize>&,
    const jaboTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline jaboTransport<Thermo, PolySize> operator*
(
    const scalar,
    const jaboTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline jaboTransport<Thermo, PolySize> operator==
(
    const jaboTransport<Thermo, PolySize>&,
    const jaboTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const jaboTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class jaboTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class jaboTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity jabo coefficients
        //  Note: input in [Pa.s], but internally uses [Pa.s/kmol]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity jabo coefficients
        //  Note: input in [W/m/K], but internally uses [W/m/K/kmol]
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline jaboTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly
        );


public:

    // Constructors

        //- Construct copy
        inline jaboTransport(const jaboTransport&);

        //- Construct as named copy
        inline jaboTransport(const word&, const jaboTransport&);

        //- Construct from Istream
        jaboTransport(Istream& is);

        //- Construct from dictionary
        jaboTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<jaboTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<jaboTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<jaboTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "jabo<" + Thermo::typeName() + '>';
        }

	inline scalar a(const int i, const scalar omega) const;

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline jaboTransport& operator=(const jaboTransport&);
        inline void operator+=(const jaboTransport&);
        inline void operator-=(const jaboTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend jaboTransport operator+ <Thermo, PolySize>
        (
            const jaboTransport&,
            const jaboTransport&
        );

        friend jaboTransport operator- <Thermo, PolySize>
        (
            const jaboTransport&,
            const jaboTransport&
        );

        friend jaboTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const jaboTransport&
        );

        friend jaboTransport operator== <Thermo, PolySize>
        (
            const jaboTransport&,
            const jaboTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const jaboTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "jaboTransportI.H"

#ifdef NoRepository
#   include "jaboTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
