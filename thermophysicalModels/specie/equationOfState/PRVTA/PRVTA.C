/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Volume-translated (Abudour's method) Peng Robinson equation of state for real gas prediction
*/

#include "PRVTA.H"
#include "IOstreams.H"
#include "Polynomial.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PRVTA<Specie>::PRVTA(Istream& is)
:
    Specie(is),
    Tc_(readScalar(is)),
    Vc_(readScalar(is)),
    Zc_(readScalar(is)),
    Pc_(readScalar(is)),
    rhoC_(readScalar(is)),
    psiSwitch_(is),
    CpCoeffs_("CpCoeffs<8>", is),
    omega_(readScalar(is)),
    rho_(0.0),
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
    is.check("PRVTA<Specie>::PRVTA(Istream& is)");
}


template<class Specie>
Foam::PRVTA<Specie>::PRVTA
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
    Vc_(readScalar(dict.subDict("equationOfState").lookup("Vc"))),
    Zc_(readScalar(dict.subDict("equationOfState").lookup("Zc"))),
    Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
    rhoC_(readScalar(dict.subDict("equationOfState").lookup("rhoC"))),
    psiSwitch_(dict.subDict("equationOfState").lookup("psiSwitch")),
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
void Foam::PRVTA<Specie>::write(Ostream& os) const
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
inline Foam::scalar Foam::PRVTA<Specie>::CpOverW(const PRVTA<Specie> pg) const
{
    return pg.CpCoeffs_/this->W();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PRVTA<Specie>& pg
)
{
    const scalar cpCoeffsOverW = CpOverW(pg);
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.Tc_
        << token::SPACE << pg.Vc_
        << token::SPACE << pg.Zc_
        << token::SPACE << pg.Pc_
	<< token::SPACE << pg.rhoC_
        << token::SPACE << pg.omega_
	<< token::SPACE << pg.rho_
	<< token::SPACE << pg.psiSwitch_
	<< token::SPACE << cpCoeffsOverW;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const PRVTA<Specie>& st)"
    );
    return os;
}

// ************************************************************************* //
