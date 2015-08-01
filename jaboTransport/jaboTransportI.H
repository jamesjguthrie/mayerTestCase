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

/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Chung's method to calculate viscosity and thermal conductivity
*/

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize>::jaboTransport
(
    const jaboTransport& pt
)
:
    Thermo(pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize>::jaboTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& muCoeffs,
    const Polynomial<PolySize>& kappaCoeffs
)
:
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize>::jaboTransport
(
    const word& name,
    const jaboTransport& pt
)
:
    Thermo(name, pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::jaboTransport<Thermo, PolySize> >
Foam::jaboTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<jaboTransport<Thermo, PolySize> >
    (
        new jaboTransport<Thermo, PolySize>(*this)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::jaboTransport<Thermo, PolySize> >
Foam::jaboTransport<Thermo, PolySize>::New(Istream& is)
{
    return autoPtr<jaboTransport<Thermo, PolySize> >
    (
        new jaboTransport<Thermo, PolySize>(is)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::jaboTransport<Thermo, PolySize> >
Foam::jaboTransport<Thermo, PolySize>::New(const dictionary& dict)
{
    return autoPtr<jaboTransport<Thermo, PolySize> >
    (
        new jaboTransport<Thermo, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::jaboTransport<Thermo, PolySize>::a
(
    const int i,
    const scalar omega
) const
{
    //only valid in nonpolar fluids
    int index = i - 1;
    scalar aZero[] = {6.32402, 0.12102e-2, 5.28346, 6.62263, 19.7454, -1.89992, 24.2745, 0.79716, -0.23816, 0.068629};
    scalar aOne[] = {50.41190, -0.11536e-2, 254.209, 38.0957, 7.63034, -12.5367, 3.44945, 1.11764, 0.067695, 0.34793};

    return aZero[index] + aOne[index] * omega;
}

template<class Thermo, int PolySize>
inline Foam::scalar Foam::jaboTransport<Thermo, PolySize>::mu
(
    const scalar p,
    const scalar T
) const
{
    //only valid in nonpolar fluids

    scalar muPolynomial = muCoeffs_.value(T)/this->W();

    scalar A = 1.16145;
    scalar B = 0.14874;
    scalar C = 0.52487;
    scalar D = 0.7732;
    scalar E = 2.16178;
    scalar F = 2.43787;
    scalar G = -6.435e-4;
    scalar H = 7.27371;
    scalar S = 18.0323;
    scalar W = -0.76830;
    scalar Tstar = (1.2593 * T)/ this->Tc_;
    scalar omegaStar = A / pow(Tstar,B) + C / exp(D*Tstar) + E / exp(F*Tstar) + G * pow(Tstar,B) * sin(S * pow(Tstar,W) - H); //reduced collision integral
    scalar omega = this->omega_;
    scalar Fc = 1 - 0.2756 * omega;
    scalar Vc = this->Vc_ / 1000; //converts to cm3/mol
    scalar sigma = 0.809 * pow(Vc,1/3);
    scalar muZero = 4.0785e-5 * (pow(this->W() * T,1/2)/(pow(Vc,2/3) * omegaStar)) * Fc;
    scalar rhoCalc = this->rho(p,T);// / this->W(); //converts to mol/litre
    scalar Y = rhoCalc * Vc / 6;
    scalar y = 0.318 * (M_PI/6) * (rhoCalc / (310));// / this->W()));
    scalar G1 = (1 - Y / 2) / pow(1-Y, 3);
    scalar expNega4Y = -a(4,omega) * Y;
    
    scalar G2 = ((a(1, omega) / y) * 
	(1 - exp(-a(4, omega) * Y))
	+ a(2,omega)*G1*exp(a(5,omega)*Y)
	+ a(3,omega)*G1)
	/ (a(1,omega)
	*a(4,omega)
	+ a(2,omega)
	+ a(3,omega));
    
    scalar muP = (36.344e-2 * pow(this->W()*this->Tc_,1/2) / pow(Vc,2/3)) * a(7,omega) * pow(Y,2) * G2 * exp(a(8,omega) + a(9,omega)/Tstar + a(10,omega)/pow(Tstar,2));
    scalar muK = muZero * (1/G2 + a(6,omega)*Y);
    scalar chungMu = muK + muP;

    //Info << "muPolynomial: " << muPolynomial << ", chungMu: " << chungMu << ", muK: " << muK << ", muP: " << muP << ", muZero: " << muZero << endl;
//    Info << Vc << endl;
    return chungMu / 10;
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::jaboTransport<Thermo, PolySize>::kappa
(
    const scalar p,
    const scalar T
) const
{
    scalar A = 1.16145;
    scalar B = 0.14874;
    scalar C = 0.52487;
    scalar D = 0.7732;
    scalar E = 2.16178;
    scalar F = 2.43787;
    scalar G = -6.435e-4;
    scalar H = 7.27371;
    scalar S = 18.0323;
    scalar W = -0.76830;
    scalar Tstar = (1.2593 * T)/ this->Tc_;
    scalar omegaStar = A / pow(Tstar,B) + C / exp(D*Tstar) + E / exp(F*Tstar) + G * pow(Tstar,B) * sin(S * pow(Tstar,W) - H); //reduced collision integral
    scalar omega = this->omega_;
    scalar Fc = 1 - 0.2756 * omega;
    scalar Vc = this->Vc_ / 1000; //converts to cm3/mol
    scalar sigma = 0.809 * pow(Vc,1/3);
    scalar muZero = 4.0785e-5 * (pow(this->W() * T,1/2)/(pow(Vc,2/3) * omegaStar)) * Fc;
    //muZero = muZero / 10;
    Info << muZero << endl;
    scalar rho = this->rho(p,T);
    scalar alpha = this->CvCalc(rho,T)/this->R()-3/2;
    scalar beta = 0.7862 - 0.7109*this->omega_ + 1.3168*pow(this->omega_,2);
    scalar Z = 2.0 + 10.5*pow((T/this->Tc_),2);
    scalar chungPsi = 1 + ((0.215+0.28288*alpha-1.061*beta+0.26665*Z)/(0.6366+beta*Z+1.061*alpha*beta));
    scalar chungKappa = 7.452*(muZero / chungPsi);
    chungKappa = chungKappa * 100;
    scalar originalKappa = kappaCoeffs_.value(T)/this->W();
    Info << chungKappa << " " << originalKappa << endl;
    return chungKappa;
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::jaboTransport<Thermo, PolySize>::alphah
(
    const scalar p, const scalar T
) const
{
    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize>&
Foam::jaboTransport<Thermo, PolySize>::operator=
(
    const jaboTransport<Thermo, PolySize>& pt
)
{
    Thermo::operator=(pt);

    muCoeffs_ = pt.muCoeffs_;
    kappaCoeffs_ = pt.kappaCoeffs_;

    return *this;
}


template<class Thermo, int PolySize>
inline void Foam::jaboTransport<Thermo, PolySize>::operator+=
(
    const jaboTransport<Thermo, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    muCoeffs_ = molr1*muCoeffs_ + molr2*pt.muCoeffs_;
    kappaCoeffs_ = molr1*kappaCoeffs_ + molr2*pt.kappaCoeffs_;
}


template<class Thermo, int PolySize>
inline void Foam::jaboTransport<Thermo, PolySize>::operator-=
(
    const jaboTransport<Thermo, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    muCoeffs_ = molr1*muCoeffs_ - molr2*pt.muCoeffs_;
    kappaCoeffs_ = molr1*kappaCoeffs_ - molr2*pt.kappaCoeffs_;
}


template<class Thermo, int PolySize>
inline void Foam::jaboTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize> Foam::operator+
(
    const jaboTransport<Thermo, PolySize>& pt1,
    const jaboTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/t.nMoles();
    scalar molr2 = pt2.nMoles()/t.nMoles();

    return jaboTransport<Thermo, PolySize>
    (
        t,
        molr1*pt1.muCoeffs_ + molr2*pt2.muCoeffs_,
        molr1*pt1.kappaCoeffs_ + molr2*pt2.kappaCoeffs_
    );
}


template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize> Foam::operator-
(
    const jaboTransport<Thermo, PolySize>& pt1,
    const jaboTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) - static_cast<const Thermo&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/t.nMoles();
    scalar molr2 = pt2.nMoles()/t.nMoles();

    return jaboTransport<Thermo, PolySize>
    (
        t,
        molr1*pt1.muCoeffs_ - molr2*pt2.muCoeffs_,
        molr1*pt1.kappaCoeffs_ - molr2*pt2.kappaCoeffs_
    );
}


template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const jaboTransport<Thermo, PolySize>& pt
)
{
    return jaboTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muCoeffs_,
        pt.kappaCoeffs_
    );
}


template<class Thermo, int PolySize>
inline Foam::jaboTransport<Thermo, PolySize> Foam::operator==
(
    const jaboTransport<Thermo, PolySize>& pt1,
    const jaboTransport<Thermo, PolySize>& pt2
)
{
    return pt2 - pt1;
}


// ************************************************************************* //
