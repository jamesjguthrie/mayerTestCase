/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014  OpenFOAM Foundation
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
    Foam::PengRobinsonGas

Description
    PengRobinsonGas gas equation of state.

SourceFiles
    PengRobinsonGasI.H
    PengRobinsonGas.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonGas_H
#define PengRobinsonGas_H

#include "autoPtr.H"
#include "Polynomial.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class PengRobinsonGas;

template<class Thermo, template<class> class Type> class thermo;

template<class Specie>
inline PengRobinsonGas<Specie> operator+
(
    const PengRobinsonGas<Specie>&,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
inline PengRobinsonGas<Specie> operator-
(
    const PengRobinsonGas<Specie>&,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
inline PengRobinsonGas<Specie> operator*
(
    const scalar,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
inline PengRobinsonGas<Specie> operator==
(
    const PengRobinsonGas<Specie>&,
    const PengRobinsonGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const PengRobinsonGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class PengRobinsonGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class PengRobinsonGas
:
    public Specie
{
    protected:

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

	//Cp calculated
	mutable scalar CpEoS_;

        //- Specific heat at constant pressure polynomial coeffs [J/(kg.K)]
        Polynomial<8> CpCoeffs_;

        //- Acentric factor [-]
        scalar omega_;

    //-Peng Robinson factors
    scalar n_; 
    scalar a0_;
    mutable scalar b_;

    mutable scalar counter;

    //CL: Variables to save the values of a, dadT and d2adT2 of the mixture
    //CL: Variables must corrected for changing temperatures
    mutable scalar aSave;
    mutable scalar daSave;
    mutable scalar d2aSave;

    //CL: save the temperature for which the save coefficients (amix,dadTmix,d2adT2mix) are correct
    mutable scalar TSave;

    //Save the Z used in the rho calculation to return the psi
    mutable scalar ZSave;

    //- Density @STD, initialise after a0, b!
    mutable scalar rhostd_;

    //Protected functions
    //CL: function updates the coefficients (aSave, daSave, d2aSave)
    inline void updateModelCoefficients(const  scalar T) const; 

public:

    // Constructors

        //- Construct from components
        inline PengRobinsonGas
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
	    const scalar& CpEoS,
	    const Polynomial<8>& CpCoeffs,
	    const scalar& omega
        );

        //- Construct from Istream
        PengRobinsonGas(Istream&);

        //- Construct from dictionary
        PengRobinsonGas(const dictionary& dict);

        //- Construct as named copy
        inline PengRobinsonGas(const word& name, const PengRobinsonGas&);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonGas> clone() const;

        // Selector from Istream
        inline static autoPtr<PengRobinsonGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<PengRobinsonGas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PengRobinsonGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

//- Return density [kg/m^3] 
        // rho0 is the starting point of the newton solver used to calculate rho
        inline scalar rho
        (
            const scalar p,
            const scalar T,
            const scalar rho0
        ) const;

	    inline scalar volTranslation(const scalar Vut, const scalar T) const;

            inline scalar rho(const scalar p,const  scalar T) const;

	    inline scalar Cp(const scalar Rho, const scalar T)const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(const scalar p, const scalar T) const;

            //- Return compression factor [-]
            inline scalar Z(const scalar p, const scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;

	    inline scalar rhostd()const;

            //CL: Model coefficient a(T)
	    inline scalar a(const  scalar T)const;

            //CL: temperature deriviative of model coefficient a(T)
	    inline scalar dadT(const  scalar T)const;

            //CL: second order temperature deriviative of model coefficient a(T)
	    inline scalar d2adT2(const  scalar T)const;

	    inline scalar a0()const;

	    inline scalar b()const;

	    inline scalar n()const;

	    inline scalar CpOverW(const PengRobinsonGas<Specie> pg) const;

	    inline scalar integral_d2pdT2_dv(const scalar rho,const scalar T) const;

	    inline scalar isobarExpCoef(const scalar rho,const scalar T)const;

	    inline scalar isothermalCompressiblity(const scalar rho,const scalar T) const;

	    inline scalar cv0(const scalar T) const;

	    inline scalar CpCalc(const scalar rho, const scalar T)const;

	    inline scalar CvCalc(const scalar rho, const scalar T)const;

            //CL: Equation of state
            inline scalar p(const scalar rho, const scalar T) const;

            //CL: first order derivatives 
            inline scalar dpdv(const scalar rho,const  scalar T) const;

            inline scalar dpdT(const scalar rho, const scalar T) const;

            inline scalar dvdT(const scalar rho,const  scalar T) const;

            inline scalar dvdp(const scalar rho, const scalar T) const;
 
	    //CL: second order derivatives, not Used At The Moment
            inline scalar d2pdv2(const scalar rho,const  scalar T) const;    		
	
            inline scalar d2pdT2(const scalar rho,const  scalar T) const;    		
          
            inline scalar d2pdvdT(const scalar rho,const  scalar T) const;    		
	
            inline scalar d2vdT2(const scalar rho,const  scalar T) const;  
	
	    inline scalar integral_p_dv(const scalar rho,const  scalar T) const;
        
        //CL: Used for Entropy
            inline scalar integral_dpdT_dv(const scalar rho,const  scalar T) const;

    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const PengRobinsonGas&);
        inline void operator-=(const PengRobinsonGas&);

        inline void operator*=(const scalar);


    // Friend operators

        friend PengRobinsonGas operator+ <Specie>
        (
            const PengRobinsonGas&,
            const PengRobinsonGas&
        );

        friend PengRobinsonGas operator- <Specie>
        (
            const PengRobinsonGas&,
            const PengRobinsonGas&
        );

        friend PengRobinsonGas operator* <Specie>
        (
            const scalar s,
            const PengRobinsonGas&
        );

        friend PengRobinsonGas operator== <Specie>
        (
            const PengRobinsonGas&,
            const PengRobinsonGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const PengRobinsonGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonGasI.H"

#ifdef NoRepository
#   include "PengRobinsonGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
