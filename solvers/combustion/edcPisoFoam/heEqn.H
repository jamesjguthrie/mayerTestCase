
    fvScalarMatrix heEqn
    (
        fvm::ddt(rho, he)
      + mvConvection->fvmDiv(phi, he)
      - fvm::laplacian(turbulence->alphaEff(), he)
      - combustion->Sh()
      - radiation->Sh(thermo) /* instead of radiation->Shs(thermo) */ 
      ==
        DpDt
      + fvOptions(rho, he)
    );

    heEqn.relax();
    fvOptions.constrain(heEqn);
    heEqn.solve();
    fvOptions.correct(he);

    //#include "limitTemperature.H"	
    
    radiation->correct();

    thermo.correct();
    
    Info<< "min/max(T) = "    << min(T).value() << ", " << max(T).value() << endl;
