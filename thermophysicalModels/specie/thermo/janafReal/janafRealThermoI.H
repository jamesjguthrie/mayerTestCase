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

#include "janafRealThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState>::janafRealThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename janafRealThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename janafRealThermo<EquationOfState>::coeffArray& lowCpCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
        lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
    }
}


template<class EquationOfState>
inline const typename Foam::janafRealThermo<EquationOfState>::coeffArray&
Foam::janafRealThermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState>::janafRealThermo
(
    const word& name,
    const janafRealThermo& jt
)
:
    EquationOfState(name, jt),
    Tlow_(jt.Tlow_),
    Thigh_(jt.Thigh_),
    Tcommon_(jt.Tcommon_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
/*        WarningIn
        (
            "janafRealThermo<EquationOfState>::limit(const scalar T) const"
        )   << "attempt to use janafRealThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;
*/
        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class EquationOfState>
inline const typename Foam::janafRealThermo<EquationOfState>::coeffArray&
Foam::janafRealThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::janafRealThermo<EquationOfState>::coeffArray&
Foam::janafRealThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    scalar cp0 = this->RR*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
    scalar cv0 = cp0 - this->RR;
    scalar rho = this->rho(p,T);
//    scalar cvReal = cv0+T*this->integral_d2pdT2_dv(rho, T);

    //scalar cpnonlimited = fabs(this->CvCalc(this->rhostd(),this->Tstd) - this->Tstd * sqr(this->dpdT(this->rhostd(), this->Tstd)) / min(this->dpdv(this->rhostd(), this->Tstd),-1));
  //  Info << "cv:  " << this->CvCalc(this->rhostd(),this->Tstd) << ", dpdT: " << this->dpdT(this->rhostd(), this->Tstd) << endl;

    //scalar cpcoeff = cpnonlimited*20; //CpCoeffs_.value(T)*20;
scalar cpReal =         std::abs
        (
            cv0
            -T*pow((this->dpdT(rho, T)),2)
                /min(this->dpdv(rho, T),-1)
        );

cpReal = 
    min
    (
        cp0*20,
        cpReal
    );
    //Info << "cp0: " << cp0 << " cpReal: " << cpReal << endl;
    return cpReal;
}

template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::cv
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T); 
    scalar cp0 = this->RR*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0]);
    scalar cv0 = cp0 - this->RR;
    scalar rho = this->rho(p,T);
    scalar cvReal = cv0+T*this->integral_d2pdT2_dv(rho, T);
    //Info << "cvReal: " << cvReal << endl; 
    return cvReal;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return this->RR*
    (
        ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
      + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    const scalar Tstd = specie::Tstd;
    return this->RR*
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janafRealThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    this->RR*
    (
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*::log(T)
      + a[6]
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::janafRealThermo<EquationOfState>::operator+=
(
    const janafRealThermo<EquationOfState>& jt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(jt);

    molr1 /= this->nMoles();
    scalar molr2 = jt.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, jt.Tlow_);
    Thigh_ = min(Thigh_, jt.Thigh_);

    if (janafRealThermo<EquationOfState>::debug && notEqual(Tcommon_, jt.Tcommon_))
    {
        FatalErrorIn
        (
            "janafRealThermo<EquationOfState>::operator+="
            "(const janafRealThermo<EquationOfState>& jt) const"
        )   << "Tcommon " << Tcommon_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << jt.Tcommon_ << " for "
            << (jt.name().size() ? jt.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<janafRealThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs_[coefLabel] =
            molr1*highCpCoeffs_[coefLabel]
          + molr2*jt.highCpCoeffs_[coefLabel];

        lowCpCoeffs_[coefLabel] =
            molr1*lowCpCoeffs_[coefLabel]
          + molr2*jt.lowCpCoeffs_[coefLabel];
    }
}


template<class EquationOfState>
inline void Foam::janafRealThermo<EquationOfState>::operator-=
(
    const janafRealThermo<EquationOfState>& jt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(jt);

    molr1 /= this->nMoles();
    scalar molr2 = jt.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, jt.Tlow_);
    Thigh_ = min(Thigh_, jt.Thigh_);

    if (janafRealThermo<EquationOfState>::debug && notEqual(Tcommon_, jt.Tcommon_))
    {
        FatalErrorIn
        (
            "janafRealThermo<EquationOfState>::operator-="
            "(const janafRealThermo<EquationOfState>& jt) const"
        )   << "Tcommon " << Tcommon_ << " for "
            << (this->name().size() ? this->name() : "others")
            << " != " << jt.Tcommon_ << " for "
            << (jt.name().size() ? jt.name() : "others")
            << exit(FatalError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<janafRealThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs_[coefLabel] =
            molr1*highCpCoeffs_[coefLabel]
          - molr2*jt.highCpCoeffs_[coefLabel];

        lowCpCoeffs_[coefLabel] =
            molr1*lowCpCoeffs_[coefLabel]
          - molr2*jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState> Foam::operator+
(
    const janafRealThermo<EquationOfState>& jt1,
    const janafRealThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;

    scalar molr1 = jt1.nMoles()/eofs.nMoles();
    scalar molr2 = jt2.nMoles()/eofs.nMoles();

    typename janafRealThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename janafRealThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        register label coefLabel=0;
        coefLabel<janafRealThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            molr1*jt1.highCpCoeffs_[coefLabel]
          + molr2*jt2.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            molr1*jt1.lowCpCoeffs_[coefLabel]
          + molr2*jt2.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janafRealThermo<EquationOfState>::debug
     && notEqual(jt1.Tcommon_, jt2.Tcommon_)
    )
    {
        FatalErrorIn
        (
            "operator+"
            "(const janafRealThermo<EquationOfState>& jt1,"
            " const janafRealThermo<EquationOfState>& jt2)"
        )   << "Tcommon " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << " != " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << exit(FatalError);
    }

    return janafRealThermo<EquationOfState>
    (
        eofs,
        max(jt1.Tlow_, jt2.Tlow_),
        min(jt1.Thigh_, jt2.Thigh_),
        jt1.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState> Foam::operator-
(
    const janafRealThermo<EquationOfState>& jt1,
    const janafRealThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs -= jt2;

    scalar molr1 = jt1.nMoles()/eofs.nMoles();
    scalar molr2 = jt2.nMoles()/eofs.nMoles();

    typename janafRealThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename janafRealThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        register label coefLabel=0;
        coefLabel<janafRealThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            molr1*jt1.highCpCoeffs_[coefLabel]
          - molr2*jt2.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            molr1*jt1.lowCpCoeffs_[coefLabel]
          - molr2*jt2.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janafRealThermo<EquationOfState>::debug
     && notEqual(jt1.Tcommon_, jt2.Tcommon_)
    )
    {
        FatalErrorIn
        (
            "operator-"
            "(const janafRealThermo<EquationOfState>& jt1,"
            " const janafRealThermo<EquationOfState>& jt2)"
        )   << "Tcommon " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << " != " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << exit(FatalError);
    }

    return janafRealThermo<EquationOfState>
    (
        eofs,
        max(jt1.Tlow_, jt2.Tlow_),
        min(jt1.Thigh_, jt2.Thigh_),
        jt1.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const janafRealThermo<EquationOfState>& jt
)
{
    return janafRealThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt),
        jt.Tlow_,
        jt.Thigh_,
        jt.Tcommon_,
        jt.highCpCoeffs_,
        jt.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::janafRealThermo<EquationOfState> Foam::operator==
(
    const janafRealThermo<EquationOfState>& jt1,
    const janafRealThermo<EquationOfState>& jt2
)
{
    return jt2 - jt1;
}


// ************************************************************************* //
