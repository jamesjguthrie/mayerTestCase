/*
    Copyright 2015, James Guthrie, University of Strathclyde
    Volume-translated (Abudour's method) Peng Robinson equation of state for real gas prediction
*/

#include "PRVTA.H"
#include "mathematicalConstants.H"
#include "thermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::PRVTA<Specie>::PRVTA
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Zc,
    const scalar& Pc,
    const scalar& rhoC,
    const string& psiSwitch,
    const Polynomial<8>& CpCoeffs,
    const scalar& omega,
    const scalar& rho
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Zc_(Zc),
    Pc_(Pc),
    rhoC_(rhoC),
    psiSwitch_(psiSwitch),
    CpCoeffs_(CpCoeffs),
    omega_(omega),
    rho_(rho),
    TSave(0),
    ZSave(0),
    aSave(0),
    daSave(0),
    d2aSave(0)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PRVTA<Specie>::PRVTA
(
    const word& name,
    const PRVTA& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Pc_(pg.Pc_),
    Vc_(pg.Vc_),
    Zc_(pg.Zc_),
    rhoC_(pg.rhoC_),
    psiSwitch_(pg.psiSwitch_),
    CpCoeffs_(pg.CpCoeffs_),
    omega_(pg.omega_),
    rho_(pg.rho_),
    a0_(pg.a0_),
    b_(pg.b_),
    n_(pg.n_),
    rhostd_(pg.rhostd_),
    TSave(0),
    ZSave(0),
    aSave(0),
    daSave(0),
    d2aSave(0)
{}


template<class Specie>
inline Foam::autoPtr<Foam::PRVTA <Specie> >
Foam::PRVTA<Specie>::clone() const
{
    return autoPtr<PRVTA<Specie> >
    (
        new PRVTA<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::PRVTA<Specie> >
Foam::PRVTA<Specie>::New
(
    Istream& is
)
{
    return autoPtr<PRVTA<Specie> >(new PRVTA<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::PRVTA<Specie> >
Foam::PRVTA<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<PRVTA<Specie> >
    (
        new PRVTA<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PRVTA<Specie>::updateModelCoefficients(const scalar T)const
{
    aSave=a0()*pow(1+n()*(1-pow(T/Tc_,0.5)),2);
    daSave=a0()*n()*(n()*sqrt(T/Tc_)-n_-1)*sqrt(T/Tc_)/T;
    d2aSave=a0()*n()*(n()+1)*sqrt(T/Tc_)/(2*pow(T,2));

    //CL: saving the temperature at which the coefficients are valid
    TSave=T;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::rhostd() const
{
    return this->rho(this->Pstd,this->Tstd,this->Pstd/(this->Tstd*this->R()));
}

//CL: Model coefficient a(T)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::a(const scalar T)const
{
    //CL: check if a has already been calculated for this temperature
    if(TSave==T)
    {
        return aSave;
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return aSave;
    }
}


//CL: temperature deriviative of model coefficient a(T)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::dadT(const scalar T)const
{
    // check if a has already been calculated for this temperature
    if(TSave==T)
    {
        return daSave;
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return daSave;
    }
}

//CL: second order temperature deriviative of model coefficient a(T)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::d2adT2(const scalar T)const
{
    // check if a has already been calculated for this temperature
    if(TSave==T)
    {
        return d2aSave;
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return d2aSave;
    }
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::a0()const
{
    return 0.457235*pow(this->RR,2)*pow(Tc_,2)/Pc_;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::b()const
{
    return 0.077796*this->RR*Tc_/Pc_;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::n()const
{
    return 0.37464+1.54226*omega_-0.26992*pow(omega_,2);
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::volTranslation(const scalar Vut, const scalar T)const
{
    scalar dpdv = (
              2*this->a(T)*
              (
                  pow(this->b(),3)-pow(this->b(),2)*Vut-this->b()*pow(Vut,2)+pow(Vut,3)
              )
              -this->RR*T*
              (
                  pow(this->b(),4)-4*pow(this->b(),3)*Vut+2*pow(this->b(),2)*pow(Vut,2)
                      +4*this->b()*pow(Vut,3)+pow(Vut,4)
              )
          )
          /(pow(this->b()-Vut,2)*pow(pow(this->b(),2)-2*this->b()*Vut-pow(Vut,2),2));

    scalar d = -pow(Vut,2)/(this->R()*Tc_) * dpdv;

    scalar deltaC = ((this->R()*Tc_)/Pc_) * (0.3074 - Zc_);

    scalar c1 = 0.4266 * Zc_ - 0.01101;

    scalar vs = (-deltaC*(0.35/(0.35+d)) + ((this->R()*Tc_)/Pc_) * (c1 - exp(-2 * d) * (c1 + 0.004)));

    return vs;

}


//returns the pressure for a given density and temperature

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::p(const scalar rho, const scalar T) const
{
    scalar Vut = this->W()/rho;
    scalar p = this->RR*T/(Vut-b())-a(T)/(pow(Vut,2)+2*b()*Vut-pow(b(),2));
    //Info << "Vut " << Vut << " Vm " << Vm << endl;
//    Info << "p " << p << endl;
    return p;
}

//Real deviative dp/dv at constant temperature 
//(molar values)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    scalar dpdv = (
              2*a(T)*
              (
                  pow(b(),3)-pow(b(),2)*Vm-b()*pow(Vm,2)+pow(Vm,3)
              )
              -this->RR*T*
              (
                  pow(b(),4)-4*pow(b(),3)*Vm+2*pow(b(),2)*pow(Vm,2)
                      +4*b()*pow(Vm,3)+pow(Vm,4)
              )
          )
          /(pow(b()-Vm,2)*pow(pow(b(),2)-2*b()*Vm-pow(Vm,2),2));

//    Info << "dpdv = " << dpdv << endl;
    return dpdv;
}

//Real deviative dp/dT at constant molar volume 
//(molar values)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::dpdT(const scalar rho, const scalar T) const
{ 
    scalar Vut = this->W()/rho;
    return this->RR/(Vut-b())-dadT(T)/(pow(Vut,2)+2*b()*Vut-pow(b(),2));
}


//Real  deviative dv/dT at constant pressure 
//by using implicit differentiation
//(molar values)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::dvdT(const scalar rho,const scalar T) const
{
    return (-1)*this->dpdT(rho, T)/this->dpdv(rho, T);
}


//(molar values)
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::dvdp(const scalar rho,const scalar T)  const
{
    scalar dvdp = 1/this->dpdv(rho, T);
    Info << "dvdp = " << dvdp << endl;
    return dvdp;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::integral_p_dv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -pow(2,0.5)*a(T)*log(this->b()*(1-pow(2,0.5))+Vm)/(4*this->b())+this->RR*T*log(Vm-this->b())+pow(2,0.5)*a(T)*log(this->b()*(pow(2,0.5)+1)+Vm)/(4*this->b());
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::Cp
(
    const scalar Rho,
    const scalar T
) const
{
    return CpCoeffs_.value(T)/this->W();
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::rho
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
 
   scalar molarVolumePrevIteration;
   scalar molarVolume;
   label iter=0;
   label maxIter_=400;
   scalar tol_=1e-8;
   scalar rho1=1500; //debug::tolerances("PRVTARhoMax", 1500); //max
   scalar rho2=1; //debug::tolerances("PRVTARhoMin", 1); //min
//if (rho0 > 105) Info << rho0<< " " << rho1 << " " << rho2 << " " << this->W() << endl;
   molarVolume=this->W()/rho0; // + volTranslation(this->W()/rho0,T);
    
//   Info << a(T) << " " << a0() << " " << b() << " " << n() << endl;
   //scalar rhoVut = p/(this->Z(p, T) * this->R() * T);
   do
   {
       molarVolumePrevIteration= molarVolume;
       label i=0;
       do
       {
       //Info << this->a0() << " " << b() << " " << n() << endl;
	//Info << molarVolume << " " << molarVolumePrevIteration << endl;
	   molarVolume=molarVolumePrevIteration-(
                   (this->p((this->W()/molarVolumePrevIteration),T) - p)
                   /(this->dpdv((this->W()/molarVolumePrevIteration),T))
               )/pow(2,i);
	//Info << molarVolume << " " << molarVolumePrevIteration << endl;
	//Info << this->p((this->W()/molarVolumePrevIteration),T) - p << endl;
	//Info << this->dpdv((this->W()/molarVolumePrevIteration),T) << endl;
           i++;
		//Info << i << endl;
           if (i>8) 
           {
               //CL: using bisection methode as backup,
               //CL: solution must be between rho=0.001 to rho=1500;
               //CL: if not, change rhoMax_ and rhoMin_
	//	Info << "if statement" << endl;
               for(i=0; i<400; i++)
               {
                   scalar f1 = this->p(rho1,T) - p;
                   scalar f2 = this->p(rho2,T) - p;
                   scalar rho3 = (rho1 + rho2)/2;
                   scalar f3 = this->p(rho3,T) - p;
	//	Info << "In the for loop " << rho0<< endl;	
                   if ((f2 < 0 && f3 > 0) || (f2 > 0 && f3 < 0))
                   {
        //		Info << "Condition 1 is true" << endl;   
	            rho1=rho3;
                   }  
                   else if ((f1 < 0 && f3 > 0)||(f1 > 0 && f3 < 0))
                   {
	//		Info << "Condition 2 is true" << endl;
                       rho2=rho3;
                   }		 
                   else	
                   {
	//		Info << "Condition 3 is true" << endl;
                       rho2=(rho2 + rho3)/2;
                   }

                   if(mag(f3) < p*tol_)
                   {
	//	       Info << "Condition 4 is true" << endl;
                       molarVolume=this->W()/rho3;
                       molarVolumePrevIteration=this->W()/rho3;
                       break; 
                   }
                   else
                   {
	//		Info << "Condition 5 is true" << endl;
                       molarVolumePrevIteration=this->W()/rho3;
                   }
	//	Info << "Still in the for loop " << rho1 << " " << rho2 << " " << rho3 << endl;
               }		  
            } 
//	Info << mag(this->p((this->W()/molarVolume), T)) << " " << mag(this->p((this->W()/molarVolumePrevIteration),T)) << endl;
	}
       while
       (
           mag(this->p((this->W()/molarVolume),T) - p)
         > mag(this->p((this->W()/molarVolumePrevIteration),T) - p)
       );

       if (iter++ > maxIter_)
       {
           FatalErrorIn
           (
               "inline scalar PRVTA::rho(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
//	Info << mag(molarVolumePrevIteration-molarVolume) << " " << tol_*(this->W()/rho0) << endl;
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

//   scalar Vt = volTranslation(molarVolume,T) + molarVolume;

//   Info << rho0/ (this->W() / molarVolume) << " " << ZSave << endl;

//    molarVolume = volTranslation(molarVolume, T) + molarVolume;

    scalar finalRho = this->W() / molarVolume; //(1 / ZSave) * (this->W() / molarVolume);

    return finalRho;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::rho(const scalar p, const scalar T) const
{
    scalar Vut = this->W()/(p/(ZSave*this->R()*T));

    scalar Vt = this->volTranslation(Vut, T) + Vut; 

    scalar rhoStart = this->W() / Vt;
    
    return this->rho(p, T, rhoStart);

/*    scalar Vut = this->W() / rhoStart;
    scalar Vt = this->volTranslation(Vut, T) + Vut;
    rhoStart = this->W() / Vt;*/

//    Info << "rhoStart: " << rhoStart << endl;
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::integral_d2pdT2_dv(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return pow(2,0.5)*d2adT2(T)*log(this->b()*(pow(2,0.5)+1)+Vm)/(4*this->b())-pow(2,0.5)*d2adT2(T)*log(this->b()*(1-pow(2,0.5))+Vm)/(4*this->b());
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::integral_dpdT_dv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -pow(2,0.5)*dadT(T)*log(this->b()*(1-pow(2,0.5))+Vm)/(4*this->b())+this->RR*log(Vm-this->b())+pow(2,0.5)*dadT(T)*log(this->b()*(pow(2,0.5)+1)+Vm)/(4*this->b());
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::isobarExpCoef(const scalar rho,const scalar T) const
{ 
    return this->dvdT(rho, T) * rho / this->W();
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::isothermalCompressiblity(const scalar rho,const scalar T) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //CL: also possible 
    //CL: return -this->dvdp(rho,T)*rho/this->W();
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::cv0(const scalar T) const
{
    scalar cp = CpCoeffs_.value(T); // /this->W();
    return cp - this->RR;    
}



template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::CvCalc(const scalar rho, const scalar T) const
{
    return this->cv0(T)+T*this->integral_d2pdT2_dv(rho, T);
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::CpCalc(const scalar rho, const scalar T) const
{
 //   scalar Vt = this->W()/rho;
    //scalar Cv = CvCalc(rho, T);
    scalar cpnonlimited = this->CvCalc(rho,T) - (T * sqr(dpdT(rho, T))) / min(dpdv(rho, T),-1);
  //  Info << "Cp: " << Cp << ", Cv: " << Cv << ", T: "<< T << endl;
return  
    min
    (
        cpnonlimited*20,
        fabs
        (
            this->CvCalc(rho,T)
            -T*pow((this->dpdT(rho, T)),2)
                /min(this->dpdv(rho, T),-1)
        )
    );
}

//drho / dp at constant H or S
template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::psi
(
    const scalar p,
    const scalar T
) const
{
    ZSave = Z(p, T);
    scalar rho = this->rho(p,T);

    scalar beta=this->isobarExpCoef(rho, T);
    scalar cp = this->CpCalc(rho,T);

    if (psiSwitch_=="psiE")
    {
    return  
    -(
        (
	    T*pow(beta,2)*(1/rho)
	    -this->isothermalCompressiblity(rho,T)*cp
	)
	/
	(
            cp*(1/rho)
            -beta*this->p(rho,T)*pow(1/rho,2)
	)
    );
    }

    else if (psiSwitch_=="psiH")
    {
    	return  
    	-(
            (T*beta*beta-beta)/cp
            -this->isothermalCompressiblity(rho,T)*rho
    	);
    }
}

template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::Z
(
    const scalar p,
    const scalar T
) const
{
    scalar a = 0.45724*sqr(this->R())*sqr(Tc_)/Pc_;
    scalar b = 0.07780*this->R()*Tc_/Pc_;
    scalar Tr = T/Tc_;

/*    scalar alpha =
	exp
	(
	    (2.0 + (0.836*Tr))
	  * (1 - pow(Tr, (0.134 + (0.508 * omega_) - (0.0467 * sqr(omega_)))))
	);
    This is the new alpha function but updateModelCoefficients need updated to comply with this */

    scalar alpha =
        sqr
        (
            1.0
          + (0.37464 + 1.54226*omega_- 0.26992*sqr(omega_))
          * (1.0 - sqrt(Tr))
        );

    scalar B = b*p/(this->R()*T);
    scalar A = a*alpha*p/sqr(this->R()*T);

    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0*B - 3.0*sqr(B);
    scalar a0 = -A*B + sqr(B) + pow3(B);

    scalar Q = (3.0*a1 - a2*a2)/9.0;
    scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54;

    scalar Q3 = Q*Q*Q;
    scalar D = Q3 + Rl*Rl;

    scalar root = -1.0;

    if (D <= 0.0)
    {
        scalar th = ::acos(Rl/sqrt(-Q3));
        scalar qm = 2.0*sqrt(-Q);
        scalar r1 = qm*cos(th/3.0) - a2/3.0;
        scalar r2 = qm*cos((th + 2.0*constant::mathematical::pi)/3.0) - a2/3.0;
        scalar r3 = qm*cos((th + 4.0*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // one root is real
        scalar D05 = sqrt(D);
	scalar S = 0;
	if ((Rl + D05) > 0) {
	    S = pow(Rl + D05, 1.0/3.0);
	}
	scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3.0;
    }

    return root;
//    return p/(this->R()*T*(p/(this->R()*T)));
}


template<class Specie>
inline Foam::scalar Foam::PRVTA<Specie>::cpMcv
(
    scalar p,
    scalar T
) const
{
    ZSave = this->Z(p,T);
  //  scalar rho = p/(ZSave*this->R()*T);
    //return this->RR * ZSave;
    return this->CpCalc(rho,T) - this->CvCalc(rho,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PRVTA<Specie>::operator+=
(
    const PRVTA<Specie>& pg
)
{
    scalar molr1 = this->nMoles();
    Specie::operator+=(pg);

    molr1 /= this->nMoles();
    scalar molr2 = pg.nMoles()/this->nMoles();

    Tc_ = molr1*Tc_ + molr2*pg.Tc_;
    Vc_ = molr1*Vc_ + molr2*pg.Vc_;
    Zc_ = molr1*Zc_ + molr2*pg.Zc_;
    Pc_ = specie::RR*Zc_*Tc_/Vc_;
    rhoC_ = molr1*rhoC_ + molr2*pg.rhoC_;
    omega_ = molr1*omega_ + molr2*pg.omega_;
    rho_ = molr1*rho_ + molr2*pg.rho_;
}


template<class Specie>
inline void Foam::PRVTA<Specie>::operator-=
(
    const PRVTA<Specie>& pg
)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(pg);

    molr1 /= this->nMoles();
    scalar molr2 = pg.nMoles()/this->nMoles();

    Tc_ = molr1*Tc_ - molr2*pg.Tc_;
    Vc_ = molr1*Vc_ - molr2*pg.Vc_;
    Zc_ = molr1*Zc_ - molr2*pg.Zc_;
    Pc_ = specie::RR*Zc_*Tc_/Vc_;
    rhoC_ = molr1*rhoC_ - molr2*pg.rhoC_;
    omega_ = molr1*omega_ - molr2*pg.omega_;
    rho_ = molr1*rho_ - molr2*pg.rho_;
}


template<class Specie>
inline void Foam::PRVTA<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PRVTA<Specie> Foam::operator+
(
    const PRVTA<Specie>& pg1,
    const PRVTA<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    const scalar Tc = molr1*pg1.Tc_ + molr2*pg2.Tc_;
    const scalar Vc = molr1*pg1.Vc_ + molr2*pg2.Vc_;
    const scalar Zc = molr1*pg1.Zc_ + molr2*pg2.Zc_;
    const scalar rhoC = molr1*pg1.rhoC_ + molr2*pg2.rhoC_;

    return PRVTA<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2),
        Tc,
        Vc,
        Zc,
        specie::RR*Zc*Tc/Vc,
	rhoC,
        molr1*pg1.omega_ + molr2*pg2.omega_,
	molr1*pg1.rho_ + molr2*pg2.rho_
    );
}


template<class Specie>
Foam::PRVTA<Specie> Foam::operator-
(
    const PRVTA<Specie>& pg1,
    const PRVTA<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    const scalar Tc = molr1*pg1.Tc_ + molr2*pg2.Tc_;
    const scalar Vc = molr1*pg1.Vc_ + molr2*pg2.Vc_;
    const scalar Zc = molr1*pg1.Zc_ + molr2*pg2.Zc_;
    const scalar rhoC = molr1*pg1.rhoC_ + molr2*pg2.rhoC_;

    return PRVTA<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2),
        Tc,
        Vc,
        Zc,
        specie::RR*Zc*Tc/Vc,
	rhoC,
        molr1*pg1.omega_ - molr2*pg2.omega_,
	molr1*pg1.rho_ - molr2*pg2.rho_
    );
}


template<class Specie>
Foam::PRVTA<Specie> Foam::operator*
(
    const scalar s,
    const PRVTA<Specie>& pg
)
{
    return PRVTA<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
	pg.rhoC_,
	pg.psiSwitch_,
	pg.CpCoeffs_,
        pg.omega_,
	pg.rho_
    );
}


template<class Specie>
Foam::PRVTA<Specie> Foam::operator==
(
    const PRVTA<Specie>& pg1,
    const PRVTA<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
