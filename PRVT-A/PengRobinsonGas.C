/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "IOstreams.H"
#include "Polynomial.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)
:
    Specie(is),
    Tc_(readScalar(is)),
    Vc_(readScalar(is)),
    Zc_(readScalar(is)),
    Pc_(readScalar(is)),
    CpEoS_(0.0),
    CpCoeffs_("CpCoeffs<8>", is),
    omega_(readScalar(is)),
    a0_(0.457235*pow(this->RR,2)*pow(Tc_,2)/Pc_),
    b_(0.077796*this->RR*Tc_/Pc_), 
    n_(0.37464+1.54226*omega_-0.26992*pow(omega_,2)),
    TSave(0.0),
    rhostd_(this->rho(this->Pstd,this->Tstd,this->Pstd/(this->Tstd*this->R()))),
    ZSave(0.0),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    counter(0.0)
{
    CpCoeffs_ *= this->W();
    is.check("PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)");
}


template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
    Vc_(readScalar(dict.subDict("equationOfState").lookup("Vc"))),
    Zc_(readScalar(dict.subDict("equationOfState").lookup("Zc"))),
    Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
    CpCoeffs_
    (
        dict.subDict("thermodynamics").lookup
        (
            "CpCoeffs<8>"
        )
    ),
    omega_(readScalar(dict.subDict("equationOfState").lookup("omega")))
    {
	CpCoeffs_ *= this->W();
    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinsonGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("thermodynamics");
    dict.add
    (
        word("CpCoeffs<8>"),
        CpCoeffs_/this->W()
    );
    os  << indent << dict.dictName() << dict;
}

template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::CpOverW(const PengRobinsonGas<Specie> pg) const
{
    return pg.CpCoeffs_/this->W();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonGas<Specie>& pg
)
{
    const scalar cpCoeffsOverW = CpOverW(pg);
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.Tc_
        << token::SPACE << pg.Vc_
        << token::SPACE << pg.Zc_
        << token::SPACE << pg.Pc_
        << token::SPACE << pg.omega_
	<< token::SPACE << pg.CpEoS_
	<< token::SPACE << cpCoeffsOverW;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const PengRobinsonGas<Specie>& st)"
    );
    return os;
}

// ************************************************************************* //
