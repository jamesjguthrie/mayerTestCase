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

/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Calculate cp, cv, etc. directly    
*/

#include "CpFromEoSThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize>::CpFromEoSThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const scalar Sf,
    const Polynomial<PolySize>& CpCoeffs,
    const typename Polynomial<PolySize>::intPolyType& hCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    EquationOfState(pt),
    Hf_(Hf),
    Sf_(Sf),
    CpCoeffs_(CpCoeffs),
    hCoeffs_(hCoeffs),
    sCoeffs_(sCoeffs)

{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize>::CpFromEoSThermo
(
    const CpFromEoSThermo& pt
)
:
    EquationOfState(pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize>::CpFromEoSThermo
(
    const word& name,
    const CpFromEoSThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_),
    Sf_(pt.Sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::e0
(
    const scalar T
) const
{
    return hCoeffs_.value(T) - this->RR*T;
}

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::cp
(
    const scalar p, const scalar T
) const
{
//    return CpCoeffs_.value(T);

    scalar rho = p/(this->Z(p,T) * this->R() * T);

    return this->CpCalc(rho,T);
}

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::cv
(
    const scalar p, const scalar T
) const
{
    scalar rho = p/(this->Z(p,T) * this->R() * T);
 //   return CpCoeffs_.value(T) - this->RR;
    return this->CvCalc(rho,T);
}

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::e
(
    const scalar rho, const scalar T
) const
{
    return (
    -this->Tstd*(this->integral_dpdT_dv(this->rhostd(),this->Tstd))
	+(this->integral_p_dv(this->rhostd(),this->Tstd))
	+this->e0(T)-this->e0(this->Tstd)
	+T*this->integral_dpdT_dv(rho,T)
	-this->integral_p_dv(rho,T));
}

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::ha
(
    const scalar p, const scalar T
) const
{
 //   scalar rhoCalc = p/(this->Z(p,T)* this->R() * T);

   // return this->e(rhoCalc,T)+p/rhoCalc*this->W()-this->Pstd/this->rhostd()*this->W();

   return hCoeffs_.value(T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::hs
(
    const scalar p, const scalar T
) const
{
    return ha(p, T) - hc();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::hc()
const
{
    return Hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::CpFromEoSThermo<EquationOfState, PolySize>::s
(
    const scalar p,
    const scalar T
) const
{
 /*   scalar rhoCalc = p/(this->Z(p,T) * this->R() * T);

    return -(this->integral_dpdT_dv(this->rhostd(),this->Tstd))
    +(this->sCoeffs_.value(T)-sCoeffs_.value(this->Tstd))
    -this->RR*log(T/this->Tstd)
    + this->integral_dpdT_dv(rhoCalc,T);
*/
   return sCoeffs_.value(T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize>&
Foam::CpFromEoSThermo<EquationOfState, PolySize>::operator=
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt
)
{
    EquationOfState::operator=(pt);

    Hf_ = pt.Hf_;
    Sf_ = pt.Sf_;
    CpCoeffs_ = pt.CpCoeffs_;
    hCoeffs_ = pt.hCoeffs_;
    sCoeffs_ = pt.sCoeffs_;
    return *this;
}


template<class EquationOfState, int PolySize>
inline void Foam::CpFromEoSThermo<EquationOfState, PolySize>::operator+=
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ + molr2*pt.Hf_;
    Sf_ = molr1*Sf_ + molr2*pt.Sf_;
    CpCoeffs_ = molr1*CpCoeffs_ + molr2*pt.CpCoeffs_;
    hCoeffs_ = molr1*hCoeffs_ + molr2*pt.hCoeffs_;
    sCoeffs_ = molr1*sCoeffs_ + molr2*pt.sCoeffs_;

}


template<class EquationOfState, int PolySize>
inline void Foam::CpFromEoSThermo<EquationOfState, PolySize>::operator-=
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ - molr2*pt.Hf_;
    Sf_ = molr1*Sf_ - molr2*pt.Sf_;
    CpCoeffs_ = molr1*CpCoeffs_ - molr2*pt.CpCoeffs_;
    hCoeffs_ = molr1*hCoeffs_ - molr2*pt.hCoeffs_;
    sCoeffs_ = molr1*sCoeffs_ - molr2*pt.sCoeffs_;
}


template<class EquationOfState, int PolySize>
inline void Foam::CpFromEoSThermo<EquationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize> Foam::operator+
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt1,
    const CpFromEoSThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;
    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();

    return CpFromEoSThermo<EquationOfState, PolySize>
    (
        eofs,
        molr1*pt1.Hf_ + molr2*pt2.Hf_,
        molr1*pt1.Sf_ + molr2*pt2.Sf_,
        molr1*pt1.CpCoeffs_ + molr2*pt2.CpCoeffs_,
        molr1*pt1.hCoeffs_ + molr2*pt2.hCoeffs_,
        molr1*pt1.sCoeffs_ + molr2*pt2.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize> Foam::operator-
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt1,
    const CpFromEoSThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs -= pt2;

    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();

    return CpFromEoSThermo<EquationOfState, PolySize>
    (
        eofs,
        molr1*pt1.Hf_ - molr2*pt2.Hf_,
        molr1*pt1.Sf_ - molr2*pt2.Sf_,
        molr1*pt1.CpCoeffs_ - molr2*pt2.CpCoeffs_,
        molr1*pt1.hCoeffs_ - molr2*pt2.hCoeffs_,
        molr1*pt1.sCoeffs_ - molr2*pt2.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const CpFromEoSThermo<EquationOfState, PolySize>& pt
)
{
    return CpFromEoSThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_,
        pt.Sf_,
        pt.CpCoeffs_,
        pt.hCoeffs_,
        pt.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::CpFromEoSThermo<EquationOfState, PolySize> Foam::operator==
(
    const CpFromEoSThermo<EquationOfState, PolySize>& pt1,
    const CpFromEoSThermo<EquationOfState, PolySize>& pt2
)
{
    return pt2 - pt1;
}


// ************************************************************************* //
