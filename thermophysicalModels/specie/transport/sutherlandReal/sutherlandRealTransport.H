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
    Foam::sutherlandRealTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    sutherlandRealTransportI.H
    sutherlandRealTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandRealTransport_H
#define sutherlandRealTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class sutherlandRealTransport;

template<class Thermo>
inline sutherlandRealTransport<Thermo> operator+
(
    const sutherlandRealTransport<Thermo>&,
    const sutherlandRealTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealTransport<Thermo> operator-
(
    const sutherlandRealTransport<Thermo>&,
    const sutherlandRealTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealTransport<Thermo> operator*
(
    const scalar,
    const sutherlandRealTransport<Thermo>&
);

template<class Thermo>
inline sutherlandRealTransport<Thermo> operator==
(
    const sutherlandRealTransport<Thermo>&,
    const sutherlandRealTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandRealTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class sutherlandRealTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sutherlandRealTransport
:
    public Thermo
{
    // Private data

        // Sutherland's coefficients
        scalar As_, Ts_;
	scalar muR_;
	scalar k_;

    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline sutherlandRealTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts,
	    const scalar muR,
	    const scalar k
        );

        //- Construct from two viscosities
        inline sutherlandRealTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2,
	    const scalar muR, const scalar k
        );

        //- Construct as named copy
        inline sutherlandRealTransport(const word&, const sutherlandRealTransport&);

        //- Construct from Istream
        sutherlandRealTransport(Istream&);

        //- Construct from dictionary
        sutherlandRealTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sutherlandRealTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<sutherlandRealTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<sutherlandRealTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherlandReal<" + Thermo::typeName() + '>';
        }

	inline scalar a(const int i, const scalar omega, const scalar muR, const scalar k) const;
        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

	inline scalar b(const int i, const scalar omega, const scalar muR, const scalar k) const;
        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline sutherlandRealTransport& operator=(const sutherlandRealTransport&);

        inline void operator+=(const sutherlandRealTransport&);

        inline void operator-=(const sutherlandRealTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend sutherlandRealTransport operator+ <Thermo>
        (
            const sutherlandRealTransport&,
            const sutherlandRealTransport&
        );

        friend sutherlandRealTransport operator- <Thermo>
        (
            const sutherlandRealTransport&,
            const sutherlandRealTransport&
        );

        friend sutherlandRealTransport operator* <Thermo>
        (
            const scalar,
            const sutherlandRealTransport&
        );

        friend sutherlandRealTransport operator== <Thermo>
        (
            const sutherlandRealTransport&,
            const sutherlandRealTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sutherlandRealTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandRealTransportI.H"

#ifdef NoRepository
#   include "sutherlandRealTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
