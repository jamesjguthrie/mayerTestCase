/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Calculate cp, cv, etc. directly    
*/

#include "CpFromEoSThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::CpFromEoSThermo<EquationOfState, PolySize>::CpFromEoSThermo
(
    Istream& is
)
:
    EquationOfState(is),
    Hf_(readScalar(is)),
    Sf_(readScalar(is)),
    CpCoeffs_("CpCoeffs<" + Foam::name(PolySize) + '>', is),
    hCoeffs_(),
    sCoeffs_(),
    CpEoS_(),
    counter(0.0)
{
    Hf_ *= this->W();
    Sf_ *= this->W();
    CpCoeffs_ *= this->W();

    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();
    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(specie::Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(specie::Tstd);
}


template<class EquationOfState, int PolySize>
Foam::CpFromEoSThermo<EquationOfState, PolySize>::CpFromEoSThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf"))),
    Sf_(readScalar(dict.subDict("thermodynamics").lookup("Sf"))),
    CpCoeffs_
    (
        dict.subDict("thermodynamics").lookup
        (
            "CpCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    hCoeffs_(),
    sCoeffs_(),
    CpEoS_()
{
    Hf_ *= this->W();
    Sf_ *= this->W();
    CpCoeffs_ *= this->W();

    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(specie::Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(specie::Tstd);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::CpFromEoSThermo<EquationOfState, PolySize>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Hf", Hf_/this->W());
    dict.add("Sf", Sf_/this->W());
    dict.add
    (
        word("CpCoeffs<" + Foam::name(PolySize) + '>'),
        CpCoeffs_/this->W()
    );
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CpFromEoSThermo<EquationOfState, PolySize>& pt
)
{
    os  << static_cast<const EquationOfState&>(pt) << tab
        << pt.Hf_/pt.W() << tab
        << pt.Sf_/pt.W() << tab
        << "CpCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.CpCoeffs_/pt.W();

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const CpFromEoSThermo<EquationOfState, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
