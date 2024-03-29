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

#include "sutherlandRealTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandRealTransport<Thermo>::sutherlandRealTransport(Istream& is)
:
    Thermo(is),
    As_(readScalar(is)),
    Ts_(readScalar(is)),
    muR_(readScalar(is)),
    k_(readScalar(is))
{
    is.check("sutherlandRealTransport<Thermo>::sutherlandRealTransport(Istream&)");
}


template<class Thermo>
Foam::sutherlandRealTransport<Thermo>::sutherlandRealTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(readScalar(dict.subDict("transport").lookup("As"))),
    Ts_(readScalar(dict.subDict("transport").lookup("Ts"))),
    muR_(readScalar(dict.subDict("transport").lookup("muR"))),
    k_(readScalar(dict.subDict("transport").lookup("k")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandRealTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("As", As_);
    dict.add("Ts", Ts_);
    dict.add("muR", muR_);
    dict.add("k", k_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandRealTransport<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st) << tab << st.As_ << tab << st.Ts_ << tab << st.muR_ << tab << st.k_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const sutherlandRealTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
