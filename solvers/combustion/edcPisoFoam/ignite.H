

// Rev. 6Dec2011

if (ign.ignited())
  {
      forAll(ign.sites(), i)
      {
          const ignitionSite& ignSite = ign.sites()[i];
          if (ignSite.igniting())
          {
	      Info << "Igniting ... " << endl;
              forAll(ignSite.cells(), icelli)
              {
                  label ignCell = ignSite.cells()[icelli];
                  if (T[ignCell] < 2000.0)
                  {
                      //Info << "Igniting cell " << ignCell << endl;
              
                      // do whatever....here we just increase the enthalpy in the cell
                      scalar magH = mag(he[ignCell]); /* scalar magH = mag(hs[ignCell]); */ 
                      he[ignCell] += ignSite.strength()*magH; /* instead of hs[ignCell] += ignSite.strength()*magH;*/ 
                  }
              }
          }
      }
  }


