
   forAll(Y, i)
    {
       volScalarField& Yi = Y[i];
       //Info<< "init Yi equation ... " << i << endl;

       fvScalarMatrix YiEqn
       (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          - fvm::laplacian(DY*rho + turbulence->mut()/Sct, Yi)
          //- fvm::laplacian(turbulence->muEff(), Yi)
          ==
            combustion->R(Yi)
          + fvOptions(rho, Yi) // added
        );

        //Info<< "solve Yi equation ... " << i << endl;	
        YiEqn.relax();
        fvOptions.constrain(YiEqn);  // added
        YiEqn.solve(mesh.solver("Yi"));
        fvOptions.correct(Yi);    // added

        Yi.max(0.0);

    }

/*

   tmp<fv::convectionScheme<scalar> > mvConvection
   (
       fv::convectionScheme<scalar>::New
       (
           mesh,
           fields,
           phi,
           mesh.divScheme("div(phi,Yi_h)")
       )
   );

   {
       reaction->correct();
       dQ = reaction->dQ();
       label inertIndex = -1;
       volScalarField Yt(0.0*Y[0]);

       forAll(Y, i)
       {
           if (Y[i].name() != inertSpecie)
           {
               volScalarField& Yi = Y[i];

               fvScalarMatrix YiEqn
               (
                   fvm::ddt(rho, Yi)
                 + mvConvection->fvmDiv(phi, Yi)
                 - fvm::laplacian(turbulence->muEff(), Yi)
                ==
                   reaction->R(Yi)
                 + fvOptions(rho, Yi)
               );

               YiEqn.relax();

               fvOptions.constrain(YiEqn);

               YiEqn.solve(mesh.solver("Yi"));

               fvOptions.correct(Yi);

               Yi.max(0.0);
               Yt += Yi;
           }
           else
           {
               inertIndex = i;
           }
       }

       Y[inertIndex] = scalar(1) - Yt;
       Y[inertIndex].max(0.0);
   }

   */
