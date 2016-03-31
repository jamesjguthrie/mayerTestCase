/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Chung's method to calculate viscosity and thermal conductivity
*/

#ifndef chungTransport_H
#define chungTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class chungTransport;

template<class Thermo, int PolySize>
inline chungTransport<Thermo, PolySize> operator+
(
    const chungTransport<Thermo, PolySize>&,
    const chungTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline chungTransport<Thermo, PolySize> operator-
(
    const chungTransport<Thermo, PolySize>&,
    const chungTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline chungTransport<Thermo, PolySize> operator*
(
    const scalar,
    const chungTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline chungTransport<Thermo, PolySize> operator==
(
    const chungTransport<Thermo, PolySize>&,
    const chungTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const chungTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class chungTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class chungTransport
:
    public Thermo
{
    // Private data

        //- Dynamic viscosity chung coefficients
        //  Note: input in [Pa.s], but internally uses [Pa.s/kmol]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity chung coefficients
        //  Note: input in [W/m/K], but internally uses [W/m/K/kmol]
        Polynomial<PolySize> kappaCoeffs_;

	scalar muR_;

	scalar k_;

    // Private Member Functions

        //- Construct from components
        inline chungTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly,
	    const scalar& muR,
	    const scalar& k
        );


public:

    // Constructors

        //- Construct copy
        inline chungTransport(const chungTransport&);

        //- Construct as named copy
        inline chungTransport(const word&, const chungTransport&);

        //- Construct from Istream
        chungTransport(Istream& is);

        //- Construct from dictionary
        chungTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<chungTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<chungTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<chungTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "chung<" + Thermo::typeName() + '>';
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

        inline chungTransport& operator=(const chungTransport&);
        inline void operator+=(const chungTransport&);
        inline void operator-=(const chungTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend chungTransport operator+ <Thermo, PolySize>
        (
            const chungTransport&,
            const chungTransport&
        );

        friend chungTransport operator- <Thermo, PolySize>
        (
            const chungTransport&,
            const chungTransport&
        );

        friend chungTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const chungTransport&
        );

        friend chungTransport operator== <Thermo, PolySize>
        (
            const chungTransport&,
            const chungTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const chungTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungTransportI.H"

#ifdef NoRepository
#   include "chungTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
