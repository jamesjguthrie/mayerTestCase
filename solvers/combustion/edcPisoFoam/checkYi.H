
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 
    // rev. 20 August 2014
    // Check if all species mass fractions are defined correctly, etc. sum Yi = 0 or sum Yi = 1 at all boundaries 
    // 
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    const fvPatchList& boundaries = mesh.boundary();
    //scalar patchArea = 0.0;

    Info << "\nCheck Mass fractions at all patches ...  " <<  endl;

    forAll(boundaries, patchI)
    { 	
        
        const fvPatch& cPatch = boundaries[patchI]; //indexed boundary definition (current patch)
     
       // if (boundaries[patchI].type() != "wall" &&  boundaries[patchI].type() != "wedge" && boundaries[patchI].type() != "processor")
        if (boundaries[patchI].type() == "inlet")
        {

	   forAll(cPatch,faceI) 
           {

   	        scalar ystar = 0.0;	   

		forAll(Y,i)
 		    ystar +=  Y[i].boundaryField()[patchI][faceI] ; 
							

	   if (ystar>0.0 && ystar<1.0)
	   {
		FatalErrorIn("edc2PisoFoam")  << "Sum of species mass fractions in " << boundaries[patchI].name()  << " is not equat to 1.0" << exit(FatalError);
	   }
					
			   
           } // end for individual boundary 

           Info << "Sum of Yi for " << boundaries[patchI].name() << " is ok ... " <<  endl;

       } // end if 

    } // end for boundaries 

