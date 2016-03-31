/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Chung's method to calculate viscosity and thermal conductivity
*/

#include "chungTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::chungTransport<Thermo, PolySize>::chungTransport(Istream& is)
:
    Thermo(is),
    muCoeffs_("muCoeffs<" + Foam::name(PolySize) + '>', is),
    kappaCoeffs_("kappaCoeffs<" + Foam::name(PolySize) + '>', is),
    muR_("muR", is),
    k_("k", is)
{
    muCoeffs_ *= this->W();
    kappaCoeffs_ *= this->W();
    muR_;
    k_;
}


template<class Thermo, int PolySize>
Foam::chungTransport<Thermo, PolySize>::chungTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    muCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "muCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    kappaCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "kappaCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    muR_
    (
	readScalar(dict.subDict("transport").lookup
	(
	    "muR"
	))
    ),
    k_
    (
	readScalar(dict.subDict("transport").lookup
	(
	    "k"
	))
    )
{
    muCoeffs_ *= this->W();
    kappaCoeffs_ *= this->W();
    muR_;
    k_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::chungTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add
    (
        word("muCoeffs<" + Foam::name(PolySize) + '>'),
        muCoeffs_/this->W()
    );
    dict.add
    (
        word("kappaCoeffs<" + Foam::name(PolySize) + '>'),
        kappaCoeffs_/this->W()
    );
    
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const chungTransport<Thermo, PolySize>& pt
)
{
    os  << static_cast<const Thermo&>(pt) << tab
        << "muCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.muCoeffs_/pt.W() << tab
        << "kappaCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.kappaCoeffs_/pt.W();

    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const chungTransport<Thermo, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
