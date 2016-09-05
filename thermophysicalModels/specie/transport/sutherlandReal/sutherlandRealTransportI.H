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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandRealTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;

    muR_ = 1.0;
    k_ = 1.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo>::sutherlandRealTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts,
    const scalar muR,
    const scalar k
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts),
    muR_(muR),
    k_(k)
{}


template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo>::sutherlandRealTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2,
    const scalar muR, const scalar k
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo>::sutherlandRealTransport
(
    const word& name,
    const sutherlandRealTransport& st
)
:
    Thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_),
    muR_(st.muR_),
    k_(st.k_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealTransport<Thermo> >
Foam::sutherlandRealTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandRealTransport<Thermo> >
    (
        new sutherlandRealTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealTransport<Thermo> >
Foam::sutherlandRealTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<sutherlandRealTransport<Thermo> >
    (
        new sutherlandRealTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandRealTransport<Thermo> >
Foam::sutherlandRealTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandRealTransport<Thermo> >
    (
        new sutherlandRealTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandRealTransport<Thermo>::a
(
    const int i,
    const scalar omega,
    const scalar muR,
    const scalar k
) const
{
    //only valid in nonpolar fluids
    int index = i - 1;
    scalar aZero[] = {6.32402, 0.12102e-2, 5.28346, 6.62263, 19.7454, -1.89992, 24.2745, 0.79716, -0.23816, 0.068629};
    scalar aOne[] = {50.41190, -0.11536e-2, 254.209, 38.0957, 7.63034, -12.5367, 3.44945, 1.11764, 0.067695, 0.34793};
    scalar aTwo[] = {-51.68010, -0.62571e-2, -168.481, -8.46414, -14.3544, 4.98529, -11.29130, 0.12348e-1, -0.81630, 0.59256};
    scalar aThree[] = {1189.02, 0.37283e-1, 3898.27, 31.4178, 31.5267, -18.1507, 69.3466, -4.11661, 4.02528, -0.72663};

    return aZero[index] + aOne[index] * omega + aTwo[index] * pow(muR, 4) + aThree[index] * k;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandRealTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
//    return As_*::sqrt(T)/(1.0 + Ts_/T);
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
    scalar Tstar = 1.2593 * (T/ this->Tc_);
    scalar omegaStar = A / pow(Tstar,B) + C / exp(D*Tstar) + E / exp(F*Tstar) + G * pow(Tstar,B) * sin(S * pow(Tstar,W) - H); //reduced collision integral
    scalar omega = this->omega_;
    scalar Fc = 1 - 0.2756 * omega + 0.059035 * pow(muR_,4) + k_;
    scalar Vc = this->Vc_ / 1000; //converts to cm3/mol
    scalar sigma = 0.809 * pow(Vc,1/3);
    scalar muZero = 4.0785e-5 * (pow(this->W() * T,1/2)/(pow(Vc,2/3) * omegaStar)) * Fc;
    scalar rhoCalc = this->rho(p,T);// / 1000; //converts to mol/litre
    scalar Y = rhoCalc * Vc / 6;
    scalar y = 0.318 * (M_PI/6) * (rhoCalc / this->rhoC_);
    scalar G1 = (1 - Y / 2) / pow(1-Y, 3);
    scalar expNega4Y = -a(4,omega, muR_, k_) * Y;

    scalar G2 = ((a(1, omega, muR_, k_) / y) *
        (1 - exp(-a(4, omega, muR_, k_) * Y))
        + a(2,omega,muR_,k_)*G1*exp(a(5,omega,muR_,k_)*Y)
        + a(3,omega,muR_,k_)*G1)
        / (a(1,omega,muR_,k_)
        *a(4,omega,muR_,k_)
        + a(2,omega,muR_,k_)
        + a(3,omega,muR_,k_));
    scalar muP = (36.344e-4 * pow(this->W()*this->Tc_,1/2) / pow(Vc,2/3)) * a(7,omega,muR_,k_) * pow(Y,2) * G2 * exp(a(8,omega,muR_,k_) + a(9,omega,muR_,k_)/Tstar + a(10,omega,muR_,k_)/pow(Tstar,2));
    scalar multiplier = 1;
    multiplier = this->rhoC_ / this->Zc_;
    //Info << multiplier << endl;
    /*if (this->W() < 75 && this->W() > 35) {
        multiplier = 1000 / (this->W() / 2);
    }
    else if (this->W() <= 35) {
        multiplier = 1000 / (this->W()-28);
    }
    else if (this->W() >= 75) {
        multiplier = 1000 / this->W();
    }*/
    if (multiplier > 100) multiplier = 100;
    //Info << multiplier << endl;
    muP *= multiplier; //100;// / (this->W());
    scalar muK = muZero * (1/G2 + a(6,omega,muR_,k_)*Y);
    scalar chungMu = (muK + muP)/10;
    return chungMu;

}

template<class Thermo>
inline Foam::scalar Foam::sutherlandRealTransport<Thermo>::b
(
    const int i,
    const scalar omega,
    const scalar muR,
    const scalar k
) const
{
    int index = i - 1;
    scalar bZero[] = {2.41657, -0.50924, 6.61069, 14.54250, 0.79274, -5.86340, 81.171};
    scalar bOne[] = {0.74824, -1.50936, 5.62073, -8.91387, 0.82019, 12.80050, 114.158};
    scalar bTwo[] = {-0.91858, -49.9912, 64.7599, -5.63794, -0.69369, 9.58926, -60.841};
    scalar bThree[] = {121.721, 69.9834, 27.0389, 74.3435, 6.31734, -65.5292, 466.775};
    return bZero[index] + bOne[index]*omega + bTwo[index] * pow(muR, 4) + bThree[index] * k;
}

template<class Thermo>
inline Foam::scalar Foam::sutherlandRealTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
 /*   scalar Cv_ = this->Cv(p, T);
    scalar kappa0 = mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
*/
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
    scalar Fc = 1 - 0.2756 * omega  + 0.059035 * pow(muR_,4) + k_;
    scalar Vc = this->Vc_ / 1000; //converts to cm3/mol
    scalar sigma = 0.809 * pow(Vc,1/3);
    scalar muZero = 4.0785e-5 * (pow(this->W() * T,1/2)/(pow(Vc,2/3) * omegaStar)) * Fc;
    //muZero *= 100;
    //Info << muZero << endl;
    scalar rhoCalc = this->rho(p,T);
    scalar Y = rhoCalc * Vc / 6;
    scalar y = 0.318 * (M_PI/6) * (rhoCalc / this->rhoC_);
    scalar G1 = (1 - Y / 2) / pow(1-Y, 3);
    scalar expNega4Y = -a(4,omega,muR_,k_) * Y;

    scalar G2 = ((a(1, omega,muR_,k_) / y) *
        (1 - exp(-a(4, omega,muR_,k_) * Y))
        + a(2,omega,muR_,k_)*G1*exp(a(5,omega,muR_,k_)*Y)
        + a(3,omega,muR_,k_)*G1)
        / (a(1,omega,muR_,k_)
        *a(4,omega,muR_,k_)
        + a(2,omega,muR_,k_)
        + a(3,omega,muR_,k_));
    scalar cv = this->cp(p,T) - this->RR;
    scalar alpha = cv/this->R()-3/2;
    scalar beta = 0.7862 - 0.7109*omega + 1.3168*pow(omega,2);
    scalar Z = 2.0 + 10.5*pow((T/this->Tc_),2);
    scalar chungPsi = 1 + alpha * ((0.215+0.28288*alpha-1.061*beta+0.26665*Z)/(0.6366+beta*Z+1.061*alpha*beta));
    scalar kappaZero = 7.452*(muZero / this->W()) * chungPsi;
    //kappaZero = kappaZero * 100;
    scalar H2 = (b(1,omega,muR_,k_)*(1-exp(-b(4,omega,muR_,k_)*Y))/Y + b(2,omega,muR_,k_)*G1 * exp(b(5,omega,muR_,k_)*Y) + b(3,omega,muR_,k_)*G1) / (b(1,omega,muR_,k_) * b(4,omega,muR_,k_) + b(2,omega,muR_,k_) + b(3,omega,muR_,k_));
    scalar kappaP = (3.039e-4*pow(this->Tc_/this->W(),1/2)/pow(this->Vc_,2/3)) * b(7,omega,muR_,k_) * pow(Y,2) * H2 * pow(T/this->Tc_,1/2);
    kappaP *= 1e3;
    scalar Pmultiplier = 1;
    scalar Kmultiplier = 1;
    if (this->W() < 75 && this->W() > 35) {
        Pmultiplier = this->W()/10;
        Kmultiplier = 1000 / (this->W() / 2);
    }
    else if (this->W() <= 35) {
        Pmultiplier = this->W();
        Kmultiplier = 100;
    }
    else if (this->W() >= 75) {
        Pmultiplier = 100/this->W();
        Kmultiplier = 1000/(this->W()/2);
    }
    //scalar multiplier = 1000/denominator;
    //if (multiplier > 100) multiplier = 100;

    kappaP *= Pmultiplier;
    //kappaP *= (rhoCalc*1000)/this->rhoC_;
    scalar kappaK = kappaZero * (1/H2 + b(6,omega,muR_,k_)*Y);
    kappaK *= Kmultiplier;
    //kappaK /= T/this->Tc_;
    scalar chungKappa = (kappaP + kappaK);
    //chungKappa /= T/this->Tc_;
//    Info << "chungKappa: " << chungKappa << ", kappaP: " << kappaP << ", kappaK: " << kappaK << ", kappaZero: " << kappaZero << endl;
    return chungKappa;

}


template<class Thermo>
inline Foam::scalar Foam::sutherlandRealTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cpv(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo>&
Foam::sutherlandRealTransport<Thermo>::operator=
(
    const sutherlandRealTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    As_ = st.As_;
    Ts_ = st.Ts_;
    muR_ = st.muR_;
    k_ = st.k_;

    return *this;
}


template<class Thermo>
inline void Foam::sutherlandRealTransport<Thermo>::operator+=
(
    const sutherlandRealTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    As_ = molr1*As_ + molr2*st.As_;
    Ts_ = molr1*Ts_ + molr2*st.Ts_;
    muR_ = molr1*muR_ + molr2*st.muR_;
    k_ = molr1*k_ + molr2*st.k_;
}


template<class Thermo>
inline void Foam::sutherlandRealTransport<Thermo>::operator-=
(
    const sutherlandRealTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    As_ = molr1*As_ - molr2*st.As_;
    Ts_ = molr1*Ts_ - molr2*st.Ts_;
    muR_ = molr1*muR_ - molr2*st.muR_;
    k_ = molr1*k_ - molr2*st.k_;
}


template<class Thermo>
inline void Foam::sutherlandRealTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo> Foam::operator+
(
    const sutherlandRealTransport<Thermo>& st1,
    const sutherlandRealTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandRealTransport<Thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_,
	molr1*st1.muR_ + molr2*st2.muR_,
	molr1*st1.k_ + molr2*st2.k_
    );
}


template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo> Foam::operator-
(
    const sutherlandRealTransport<Thermo>& st1,
    const sutherlandRealTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandRealTransport<Thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_,
	molr1*st1.muR_ - molr2*st2.muR_,
	molr1*st1.k_ - molr2*st2.k_
    );
}


template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandRealTransport<Thermo>& st
)
{
    return sutherlandRealTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_,
	st.muR_,
	st.k_
    );
}


template<class Thermo>
inline Foam::sutherlandRealTransport<Thermo> Foam::operator==
(
    const sutherlandRealTransport<Thermo>& st1,
    const sutherlandRealTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
