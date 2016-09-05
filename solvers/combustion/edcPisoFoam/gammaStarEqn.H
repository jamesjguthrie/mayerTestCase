  
    if (edcLES)
    {	 

    	
	
	kSgs.internalField() =  turbulence->k();
 	eSgs.internalField() =  turbulence->epsilon() ;
	bound ( eSgs, eSgsMin_ );   
	
	

	ReTauLES_.internalField() = kSgs.internalField() * kSgs.internalField() / eSgs.internalField() / turbulence->mu() * thermo.rho() ; 
	bound ( ReTauLES_, 0 );   

	gammaL = pow(  turbulence->mu() / thermo.rho() * eSgs / pow( kSgs, 2.0)  , (0.25) );		
	tauK = pow( turbulence->mu() / thermo.rho() / eSgs , (0.5) );  // [sec]
	
	
	
	tauStar = tauK; // tauDelta ;
	bound ( tauStar, tauMin_ );   
	

   }
    
