/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "chemistryModel.H"
#include "reactingMixture.H"
//#include "UniformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::chemistryModel
(
    const fvMesh& mesh
)
:
    CompType(mesh),
    ODE(),

    // initialize RADAU5 parameters:
    RADAU5s(
        this->thermo().composition().Y().size()+1,      //number of equations to solve (=number of species + 1)
        _EDC_solverX_,            //solution vector (double array of dimension nSpecie_+1)
        0.,                 //startTime (will be reset later)
        10.,                //endTime (will be reset later)
        1.,                 //time step for intermediate output
        0,                   //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
        new scalar(readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("combustionProperties")).subDict("edcPSRCoeffs").lookup("relativeTolerance"))),
                            //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
        new scalar(readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("combustionProperties")).subDict("edcPSRCoeffs").lookup("absoluteTolerance"))),
                            //absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
        false,               //provide data for output functions
        1.e-3,              //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
        0.,                 //maximal step size ( if 0. code sets (xend - x) )
        readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("combustionProperties")).subDict("edcPSRCoeffs").lookup("maxIterations")) ,
                            //max. number of steps ( if 0. then 1e5 is assumed)
        0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
        0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
        0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed),
        0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
        false,              //Switch for the computation of the Jacobian (0= Jacobian calculated internally)
        this->thermo().composition().Y().size()+1, //Switch for the banded structure of the Jacobian
        this->thermo().composition().Y().size()+1, //Switch for the banded structure of the Jacobian
        0,                  //Gives information on the mass-matrix ( 0 = mass-matrix is unity)
        0.,                 //Switch for the banded structure of the mass-matrix
        0.,                 //Switch for the banded structure of the mass-matrix
        0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
        true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
        0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
        0,                  //Dimension of the index 2 variables. Default nind2 = 0.
        0,                  //Dimension of the index 2 variables. Default nind2 = 0.
        1,              //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default];
                            // If npred = 2  classical step size control (less safe; slightly faster)
        0,                  //Default m1 = 0
        0,                  //Default m2 = m1
        false,              //If != 0 Jacobian matrix transformed to Hessenberg form
                            //(advantageous for large systems with full Jacobian)
        0,                  //Stopping criterion for Newton's method, usually chosen < 1.
                            //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
        0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
        0.,                 //This saves, together with a large thet, lu-decompositions and
                            //computing time for large systems. for small systems one may have
                            //quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
                            //quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

        0.),                //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly


    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    _EDC_solverX_(new double[nSpecie_+1]),
    //_EDC_vfYStar_(nSpecie_),
    //_EDC_vfYStar_(this->thermo().composition().Y()),
    _EDC_sfYStar_(nSpecie_),
    _EDC_mDotStar_(mesh.nCells(),0.),
    _EDC_gammaStar_(mesh.nCells(),0.),
    _EDC_ReTau_(mesh.nCells(),0.),
   // _EDC_tau_(mesh.nCells(),0.),
    _EDC_edcFactor_(mesh.nCells(),0.),
    _EDC_RR_(nSpecie_),
    RR_(nSpecie_)
{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "RR." + Y_[fieldI].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );

    }


    forAll(_EDC_RR_, specieI)
    {
       _EDC_RR_.set(specieI, new scalarField(mesh.nCells(), 0.));
        // allocate memory for fine structure data and initialize with cell average mass fractions
       _EDC_sfYStar_.set(specieI, new scalarField( Y_[specieI].internalField() ));
        //Info << specieI << endl;

    }





    //Info<< "chemistryModel: Number of species = " << nSpecie_    << " and reactions = " << nReaction_ << endl;



    // get access to combustionPropertiesDict to read EDC parameters
    const IOdictionary& combustionProperties = this->db().objectRegistry::lookupObject<IOdictionary>("combustionProperties");

    // LES or RANS version
    // 0 - RANS
    // 1 -  LES
    _EDCtype_ = 	readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("EDCSimulationType"));

    // URANS EDC PSR constants
    _EDC_CD1_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("CD1"));
    _EDC_CD2_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("CD2"));
    _EDC_CD3_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("CD3"));

    _EDC_gammaStarClipFactor_ = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("gammaStarClipFactor"));

    //_EDC_reactorType_ = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("reactorType"));

    // read Local Extinction polynomial coefficients from chemistryProperties
    _EDC_relaxFineStructures_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("underRelaxFineStructures"));


    // LES EDC PSR constants



    // EDC temperature limits
    _EDC_TMin_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("TMin"));
    _EDC_TMax_=readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("TMax"));

    // RADAU5 properties
    // RADAU5 integration absolute tolerance
    _EDC_absoluteTolerance_   = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("absoluteTolerance"));
    // RADAU5 integration relative tolerance
    _EDC_relativeTolerance_ = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("relativeTolerance"));
    // RADAU5 integration max. allowed number of iterations
    _EDC_maxIterations_	 = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("maxIterations"));


     // use binaryTable to store and retrieve ODE integration results?
    _EDC_useBinaryTree_  = readBool(combustionProperties.subDict("edcPSRCoeffs").lookup("useBinaryTree"));

    // tolerance for retrieving values and growing points
    _EDC_tableErr_   = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("binaryTreeTolerance"));

    // tolerance for retrieving values and growing points (read as scalar to allew e.g. "1.e5")
    _EDC_tableSize_  = readScalar(combustionProperties.subDict("edcPSRCoeffs").lookup("binaryTreeSize"));

    //set these values to <0 to indicate that they have not yet been initialized
    _EDC_hsMax_ = -1;
    _EDC_mDotStarMax_ = -1;
    _EDC_gammaStarMax_ = -1;

    Info << "edcPSR/Number of species: " << nSpecie_  << " and reactions = " << nReaction_ << endl;




    // construct table with fundamental elements and their number of atoms for each specie
    EDC_constructSpecieAtomsTable();




    // initialize binaryTree (if wanted)
    if(_EDC_useBinaryTree_)
    {
        _EDC_resultTable_ = new binaryTree(label(_EDC_tableSize_));
        Info << "maxTableSize: " << label(_EDC_tableSize_) << endl;
    }




}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryModel<CompType, ThermoType>::~chemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




template<class CompType, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::chemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<scalarField> tom(new scalarField(nSpecie_, 0.0));
    scalarField& om = tom();

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return tom;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{

    const Reaction<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(nSpecie_, 0.0);
    for (label i = 0; i < nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    const scalar kf = R.kf(p, T, c2);
    const scalar kr = R.kr(kf, p, T, c2);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    Info << endl << endl << endl << "************** DERIVATIVES FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;

}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{

    Info << endl << endl << endl << "************** JACOBIAN FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;

}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::tc() const
{


    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );


    Info << endl << endl << endl << "************** TC FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
    return ttc;

}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );


    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                const scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModel<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::chemistryModel<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    //return nSpecie_ + 2;

    Info << endl << endl << endl << "************** nEQNS FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
    return 0;

}


template<class CompType, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::chemistryModel<CompType, ThermoType>::calculateRR
(
    const label reactionI,
    const label specieI
) const
{


    tmp<DimensionedField<scalar, volMesh> > tRR
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );


    Info << endl << endl << endl << "************** nEQNS FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
    return tRR;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::calculate()
{


       // CompType::correct();

        if(!this->chemistry_)
        {
            Info << "Chemistry disabled in const/chemistryProperties!!" << endl;
            return;
        }

        if 	(_EDCtype_ == 0)
        {
            // RAS modeling // Info << "EDC RAS solving ... " << endl;
            EDC_RAS_solve();
        }
        else if (_EDCtype_ == 1)
        {
            //Info << "EDC LES solving ... " << endl;
            EDC_LES_solve();
        }


        else if (_EDCtype_ == 2)
        {
            // SAS modeling
            // EDC_RAS_solve();

            Info << "EDC for SAS is not implemnted !!!" << endl;
            return;
        }



        EDC_writeCombustionTerms();




}



// ************************************************************************* //

template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::RADAU5derivative(double x, double *y, double *f) const
{
        // in case of PSR we do not use the value of x explicitly
    //Info << "edcPSR/RADAU5Derivative() ... "<< endl;

    Foam::scalarField vectorY(nSpecie_+1, 0.0);
    Foam::scalarField vectordYdt(nSpecie_+1, 0.0);

    // assign each value of array type "double" to a new array of type "scalar"
    forAll(vectorY,i)
    {
        vectorY[i] = y[i];
    }

    //Info << "edcPSR/EDC_derivative() ... "<< endl;

    this->EDC_derivative(vectorY,vectordYdt);



    // and assign back from "scalar" to "double"
    forAll(vectorY,i)
    {
        f[i]=vectordYdt[i];
    }



}


// Jacobian left empty
template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::RADAU5jacobian(double x, double *y, double **J) const
{
    Info << endl << endl << endl << "************** JACOBIAN FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
}

// Mass matrix is unity, therefore left empty
template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::RADAU5mass(double **M) const
{
    Info << endl << endl << endl << "************** MASS FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;

}



template<class CompType, class ThermoType>
void  Foam::chemistryModel<CompType, ThermoType>::EDC_constructSpecieAtomsTable()
{

   std::map<char,unsigned int> mapTmp_;


   _EDC_specieAtomsTable_.clear();


   _EDC_element_ el;
   std::vector<_EDC_element_> el_data;
   el_data.clear();


   int countC = 0;
   int countN = 0;
   int countH = 0;
   int countO = 0;
   int countAR = 0;



    _EDC_mixtureWT_ = 0.0;


    const basicMultiComponentMixture& mixture = dynamic_cast<const basicMultiComponentMixture&>(this->thermo());


    for (int r = 0; r< nSpecie_; r++)
    {
        _EDC_mixtureWT_ += specieThermo_[r].W();
    }


    for (int i = 0; i< nSpecie_; i++)
    {

    el_data.clear();
    //std::string specieName = specieThermo_[i].name();

    std::string specieName = mixture.species()[i];

    std::string::iterator it = specieName.begin();



    //el.symbol = specieName[0];
    el.symbol = *it;
    el.flag = false;
    el.index.clear();
    el.count = 1;


        mapTmp_.clear();


        while( it != specieName.end() )
    {


        if (*it == 'A') // for AR
        {

            el.symbol = *it;
            el.flag = false;
            el.index.clear();
            el.count = 1;
                    el.index += '1';
            el_data.push_back(el);
            it = specieName.end();

        }
        else if (*it == '1' || *it == '2' || *it == '3' ||  *it == '4' ||  *it == '5' ||  *it == '6' ||  *it == '7' ||  *it == '8' ||  *it == '9' || *it == '0')
        {
            el.count = 0;
            el.index += *it;

            specieName.erase(it);
            it = specieName.begin();

            if (it == specieName.end())
            {
                // if the end of formula
                el_data.push_back(el);
            }
            else if ( *it == 'C' || *it == 'H' || *it == 'N' || *it == 'O' || *it == '(' || *it == ')' || *it == 'S')
            {
                // or the next element
                el_data.push_back(el);
            }



        }
        else if (*it == '(' || *it == ')' || *it == 'S')
        {
            // do nothing for (S)
            specieName.erase(it);
            it = specieName.begin();
        }
        else
        {

            el.symbol = *it;
            el.flag = false;
            el.index.clear();
            el.count = 1;

            specieName.erase(it);
            it = specieName.begin();

            if (it == specieName.end() || ( *it == 'C' || *it == 'H' || *it == 'N' || *it == 'O' || *it == '(' || *it == ')' || *it == 'S' )  )
            {
                el.index += '1';
                el_data.push_back(el);
            }
        }



    }
     // end while




    //Info << "specie: " << specieThermo_[i].name() << "\t" << "WT: " << "\t" << specieThermo_[i].W() << endl;
    Info << "specie: " << mixture.species()[i] << "\t" << "WT: " << "\t" << specieThermo_[i].W() << endl;

    for (unsigned int j = 0; j< el_data.size(); j++)
    {
        Info << el_data[j].symbol <<  "\t" << atoi(el_data[j].index.c_str()) <<   endl;

        if  (el_data[j].symbol == 'C')
        {
            countC +=  atoi(el_data[j].index.c_str());
            mapTmp_.insert(std::pair<char,unsigned int>('C', atoi(el_data[j].index.c_str()) ) );
        }
        else if  (el_data[j].symbol == 'H')
        {
            countH +=  atoi(el_data[j].index.c_str());
            mapTmp_.insert(std::pair<char,unsigned int>('H', atoi(el_data[j].index.c_str()) ) );

        }
        else if  (el_data[j].symbol == 'O')
        {
            countO +=  atoi(el_data[j].index.c_str());
            mapTmp_.insert(std::pair<char,unsigned int>('O', atoi(el_data[j].index.c_str()) ) );

        }
        else if  (el_data[j].symbol == 'N')
        {
            countN +=  atoi(el_data[j].index.c_str());
            mapTmp_.insert(std::pair<char,unsigned int>('N', atoi(el_data[j].index.c_str()) ) );

        }
        else if  (el_data[j].symbol == 'A')
        {
            countAR +=  atoi(el_data[j].index.c_str());
            mapTmp_.insert(std::pair<char,unsigned int>('AR', atoi(el_data[j].index.c_str()) ) );

        }

    }

        _EDC_specieAtomsTable_.push_back(mapTmp_);

    } // end for specie


    const scalar WTC = 12.01115;
    const scalar WTH = 1.00797;

    const scalar WTN = 14.00674;
    const scalar WTO = 15.9994;
    const scalar WTAR = 39.948;

    const scalar mWC_ = countC * WTC;
    const scalar mWH_ = countH * WTH;
    const scalar mWN_ = countN * WTN;
    const scalar mWO_ = countO * WTO;
    const scalar mWAR_ = countAR * WTAR;

    const scalar tmpWT_ =  (mWC_ +  mWH_ + mWN_ + mWO_ + mWAR_ );

    Info << "C:" 	  << "\t" <<  mWC_ <<  endl;
    Info << "H:" 	  << "\t" <<  mWH_ <<  endl;
    Info << "N:" 	  << "\t" <<  mWN_ <<  endl;
    Info << "O:" 	  << "\t" <<  mWO_ <<  endl;
    Info << "AR:" 	  << "\t" <<  mWAR_ <<  endl;

    Info << endl << "Species information: " << endl;
    Info << "summ elem: " << "\t" <<  round(tmpWT_) <<  endl;
    Info << "mixtureWT: " << "\t" <<  round(_EDC_mixtureWT_)    <<  endl;

    if (_EDC_specieAtomsTable_.size()!= nSpecie_ )
    {
    FatalErrorIn
        (
            "combustion::chemistryModel::EDC_constructSpecieAtomsTable"
            "("")"
    )   << "Size of specieAtomsTable is not equal to N species"  << abort(FatalError);
    }
    if (round(tmpWT_) != round(_EDC_mixtureWT_) )
    {
    FatalErrorIn
        (
            "combustion::chemistryModel::EDC_constructSpecieAtomsTable"
            "("")"
    )   << "calculated total atomic mass of a mxiture  is not equal to the those form a kinetic mechanism "  << abort(FatalError);
    }


}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_RAS_solve()
{

    const scalar Cmu  = 0.09;
    const scalar X    = 1.00;

    Info << "EDC/Perfectly Stirred Reactor for RAS" << endl;



    // get flow parameters: rho, k, epsilon, nu
    const scalarField rho = this->thermo().rho();

    const scalarField& k = this->db().objectRegistry::lookupObject<volScalarField>("k").internalField();

    scalarField epsilon(rho.size(),0.);

    if (this->db().objectRegistry::foundObject<volScalarField>("epsilon"))
    {
        const scalarField& epsilonConst = this->db().objectRegistry::lookupObject<volScalarField>("epsilon").internalField();
        epsilon = epsilonConst;
    }
    // or maybe we can get omega
    else if (this->db().objectRegistry::foundObject<volScalarField>("omega"))
    {
        const scalarField& omega = this->db().objectRegistry::lookupObject<volScalarField>("omega").internalField();
        epsilon = Cmu*omega*k;
    }
    // or we are lost!!
    else
    {
        FatalError << "Could not access epsilon or omega data from turbulence model!!" << exit(FatalError);
    }
    const scalarField nu = this->thermo().mu().internalField()/rho;
    // calculate EDC model variables: gammaStar and mDotStar
    
    scalarField gammaL= _EDC_CD1_ * pow(  nu*epsilon / pow(k,2)  , ( 0.25 ) );
    _EDC_gammaStar_= pow(gammaL,2.);
    _EDC_gammaStar_= min (_EDC_gammaStarClipFactor_, _EDC_gammaStar_);

    _EDC_mDotStar_ = _EDC_CD2_ * pow( epsilon/nu , 0.5 );

    _EDC_tau_  = _EDC_CD3_ * _EDC_mDotStar_;

    // calculate reaction rate coefficient
    _EDC_edcFactor_ = rho * _EDC_gammaStar_ * _EDC_mDotStar_ * X  / (1.0 - _EDC_gammaStar_*X);

    _EDC_ReTau_ = k * k / nu / epsilon;


    // calculate fine structure composition according to fast chemistry approach
    Info << "Integrating PSR equations from 0 to steady state ... ";
    EDC_updateYStar();
    Info << "done" << endl;

    //transform reaction rate to un-normalized specie space

    /*
    volScalarField _edcFactor_
    (
            IOobject
            (
                "_edcFactor_",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
    );

    _edcFactor_.internalField() = _EDC_edcFactor_;
    */

    /*
    // OF 2.1.1 implementation
    forAll(Y_,specieI)
       RR_[specieI] = (-1.0) * _edcFactor_* (Y_[specieI] - _EDC_vfYStar_[specieI]);
    */


    forAll(rho, celli)
    {        

        for (label i=0; i<nSpecie_; i++)
        {
            //RR_[i][celli] = (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_vfYStar_[i][celli]);
                  RR_[i][celli] =  (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_sfYStar_[i][celli]);
             _EDC_RR_[i][celli] =  (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_sfYStar_[i][celli]);


              //Info << "specie\t" << i << "\tRR:" << RR_[i][celli]*1e+10 << "\tEDCRR:" << _EDC_RR_[i][celli]*1e+10  <<   endl;

        }
    }

    /*
    forAll(_EDC_RR_,specieI)
    {
       _EDC_RR_[specieI]= (-1.0) * _EDC_edcFactor_ * (Y_[specieI].internalField() - _EDC_vfYStar_[specieI].internalField());
    }
        */


    //forAll(RR_,specieI)
    //   RR_[specieI]  = _EDC_RR_[specieI];



}



template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_updateYStar()
{

    scalar startTime = 0.0;
    scalar stopTime = 10.0;

    if(_EDC_useBinaryTree_)
    {


        //USE BINARY TABULATION

                //initialize some statistical data
                label added=0;
                label grown=0;
                label retrieved=0;

                // if not yet happened, initialize normalization values for hs, mDotStar and gammaStar
                // (NB this cannot be done in constructor, since there is no turbulence yet defined [requires thermo first])

                if(_EDC_mDotStarMax_ < 0)
                {
                    //get max order of magnitude of some fields
                    _EDC_hsMax_ = Foam::max(mag(this->thermo().he())).value();
                    _EDC_mDotStarMax_ = Foam::max(_EDC_mDotStar_);
                    _EDC_gammaStarMax_ = Foam::max(_EDC_gammaStar_);
                }

                forAll(Y_[0].internalField(), cellI)
                {
                    // initialize query vector (Yi , hStar/hBar , mDotStar, gammStar)
                    scalarField queryVector(Y_.size()+3);
                    //load speciesData
                    forAll(Y_,specieI)
                    {
                        queryVector[specieI]=_EDC_sfYStar_[specieI][cellI];
                    }

                    // Info << " max Values" << hsMax_ << "   " << mDotStarMax_ << "   " << gammaStarMax_ << endl;

                    // add further values and normalize them to more or less 0..1
                    queryVector[Y_.size()] = this->thermo().he().internalField()[cellI]/_EDC_hsMax_;
                    queryVector[Y_.size()+1] = _EDC_mDotStar_[cellI]/_EDC_mDotStarMax_;
                    queryVector[Y_.size()+2] = _EDC_gammaStar_[cellI]/_EDC_gammaStarMax_;

                    //query reactionTable
                    chemPoint* tableEntry = _EDC_resultTable_->findClosest(queryVector);

                    //found something??
                    if (tableEntry != NULL && ( tableEntry->inEOA(queryVector) ))
                    {

                        scalarField resultVector = tableEntry->r();
                        //unpack result vector
                        //unpack species (and underrelax)
                        forAll(_EDC_sfYStar_,specieI)
                        {
                            _EDC_sfYStar_[specieI][cellI] = _EDC_sfYStar_[specieI][cellI]*(1.-_EDC_relaxFineStructures_) + resultVector[specieI]*_EDC_relaxFineStructures_;

                        }

                        // some statistical info
                        retrieved++;
                        //and continue with next cell
                        continue;
                    }
                    else
                    {
                        // retrieve unsuccessfull, do direct integration with and add result to table (or grow it)
                        //initialize solution vector
                        forAll(Y_,specieI)
                        {
                            _EDC_solverX_[specieI] = _EDC_sfYStar_[specieI][cellI];
                        }
                        _EDC_solverX_[nSpecie_]=cellI;

                        // set stopTime for this cell
                        stopTime=min((stopTime-startTime),100*1/(max(_EDC_mDotStar_[cellI],SMALL)));

                        // for this, reset the RADAU5 parameters
                        double atol = _EDC_absoluteTolerance_;
                        double rtol = _EDC_relativeTolerance_;

                        // tell RADAU5 the start and endTimes
                        ResetT(
                            _EDC_solverX_,            //solution vector (double array of dimension nSpecie_+1)
                            startTime,                 //startTime (will be reset later)
                            stopTime,                //endTime (will be reset later)
                            1.,                 //time step for intermediate output
                            0,                  //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
                            &atol,              //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
                            &rtol,              //absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
                            0.,                 //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
                            0.,                 //maximal step size ( if 0. code sets (xend - x) )
                            _EDC_maxIterations_ ,    //max. number of steps ( if 0. then 1e5 is assumed)
                            0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
                            0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
                            0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed),
                            0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
                            nSpecie_+1,         //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
                            nSpecie_+1,         //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
                            0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
                            0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
                            0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
                            true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
                            0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
                            0,                  //Dimension of the index 2 variables. Default nind2 = 0.
                            0,                  //Dimension of the index 2 variables. Default nind2 = 0.
                            1,                  //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default];
                                                //                            If npred = 2  classical step size control (less safe; slightly faster)
                            0,                  //Default m1 = 0
                            0,                  //Default m2 = m1
                            false,              //If != 0 Jacobian matrix transformed to Hessenberg form
                                                //(advantageous for large systems with full Jacobian)
                            0,                  //Stopping criterion for Newton's method, usually chosen < 1.
                                                //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
                            0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
                            0.,                 //This saves, together with a large thet, lu-decompositions and
                                                //computing time for large systems. for small systems one may have
                                                //quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
                                                //quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

                            0.),                //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly
                        // and solve ODEs for this cell

                        Integrate();
                        // => solverX now contains result of ODE integration

                        // do explicit underrelaxation of fine structures
                        forAll(_EDC_sfYStar_,specieI)
                        {
                            _EDC_sfYStar_[specieI][cellI]=_EDC_sfYStar_[specieI][cellI]*(1.-_EDC_relaxFineStructures_) + _EDC_solverX_[specieI]*_EDC_relaxFineStructures_;
                        }

                        //store result in table
                        scalarField resultVector(Y_.size());
                        forAll(Y_,specieI)
                        {
                            resultVector[specieI]=_EDC_sfYStar_[specieI][cellI];
                        }

                        scalarField  queryTolerances(queryVector.size(),1.);
                        scalarField  resultTolerances(resultVector.size(),1.);

                        // check if grow is possible
                        if (tableEntry!=NULL)
                        {
                            if (!tableEntry->checkSolution(queryVector,resultVector))
                            {
                                //not possible, add new leaf
                                chemPoint newEntry(queryVector,resultVector,queryTolerances,resultTolerances,_EDC_tableErr_);
                                _EDC_resultTable_->insert(newEntry);
                                added ++;
                            }
                            else
                            {
                                grown++;
                            }
                        }
                        else
                        {
                            chemPoint newEntry(queryVector,resultVector,queryTolerances,resultTolerances,_EDC_tableErr_);
                            _EDC_resultTable_->insert(newEntry);
                            added ++;
                        }
                    }   // end of integration

                }   // end of cell loop

                //and output some statistics
                Info << "table size: " << _EDC_resultTable_->size() << endl;
                Info << "added points: " << added << endl;
                Info << "retrieved points: " << retrieved << endl;
                Info << "grown points: " << grown << endl;




    }
    else
    {
        //DO NOT USE TABULATION, INTEGRATE FOR EVERY CELL
        forAll(Y_[0].internalField(), cellI)
        {
            //initialize solution vector
            forAll(Y_,specieI)
            {
                //_EDC_solverX_[specieI] = YStar_[specieI][cellI];  // old one

                _EDC_solverX_[specieI] = _EDC_sfYStar_[specieI][cellI];

                //Info << specieI << endl;
                //Info <<  specieI << "   " << solverX[specieI] << endl;
                //EDC_writeScalarField(solverX[specieI],"edc_nu");
            }

            _EDC_solverX_[nSpecie_]=cellI;
           // Info <<  nSpecie_ << "   " << solverX[nSpecie_] << endl;



            // set stopTime for this cell
            stopTime = min((stopTime-startTime),100*1/(max(_EDC_mDotStar_[cellI],SMALL)));

            // for this, reset the RADAU5 parameters
            double atol = _EDC_absoluteTolerance_;
            double rtol = _EDC_relativeTolerance_;

            // Info << "edcPSR/Reset() ... "<< endl;
            // tell RADAU5 the start and endTimes

            ResetT(
                _EDC_solverX_,      //solution vector (double array of dimension nSpecie_+1)
                0.0,          //startTime (will be reset later)
                10.0,          //endTime (will be reset later)
                1.,           //time step for intermediate output
                0,            //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
                // Always should be zero !!!!
                &atol,       //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
                &rtol, 	//absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
                0.,         //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
                0.,                 //maximal step size ( if 0. code sets (xend - x) )
                 _EDC_maxIterations_ ,  //max. number of steps ( if 0. then 1e5 is assumed)
                0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
                0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
                0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed),
                0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
                nSpecie_+1, //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
                nSpecie_+1,  //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
                0.,   //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
                0.,  //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
                0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
                true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
                0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
                0,                  //Dimension of the index 2 variables. Default nind2 = 0.
                0,                  //Dimension of the index 2 variables. Default nind2 = 0.
                1,   //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default];
                                    //   If npred = 2  classical step size control (less safe; slightly faster)
                0,                  //Default m1 = 0
                0,                  //Default m2 = m1
                false,              //If != 0 Jacobian matrix transformed to Hessenberg form
                                    //(advantageous for large systems with full Jacobian)
                0,    //Stopping criterion for Newton's method, usually chosen < 1.
                       //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
                0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
                0.,                 //This saves, together with a large thet, lu-decompositions and
                                    //computing time for large systems. for small systems one may have
                                    //quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
                                    //quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

                0.);  //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly


            // and solve ODEs for this cell
            //Info << "edcPSR/Integrate() ... "<< endl;
            Integrate();

            forAll(_EDC_sfYStar_,specieI)
            {
                //YStar_[specieI][cellI]=YStar_[specieI][cellI]*(1.-relaxFineStructures_) + solverX[specieI]*relaxFineStructures_;

                _EDC_sfYStar_[specieI][cellI] = _EDC_sfYStar_[specieI][cellI]*(1.-_EDC_relaxFineStructures_) + _EDC_solverX_[specieI]*_EDC_relaxFineStructures_;

                //_EDC_vfYStar_[specieI].internalField()[cellI] =
                //        _EDC_vfYStar_[specieI].internalField()[cellI]*(1.-_EDC_relaxFineStructures_) + _EDC_solverX_[specieI]*_EDC_relaxFineStructures_;

            }


                //Info << "edcPSR/updateYStar() ...f... "<< endl;
            }


    }


} // end EDC_updateYStar();


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_iterateTAndRho
(Foam::scalar& TResult, Foam::scalar& rho, Foam::scalar h, Foam::scalar p, Foam::scalarField& Y, Foam::scalar TInit) const
{


    ThermoType mixture(0.0*specieThermo_[0]);
    // use size from global variable, so that YStarCell (which contains the cell index as Y_.size+1)
    //forAll(Y_,specieI)
    //{
    //    mixture += Y[specieI]/specieThermo_[specieI].W()*specieThermo_[specieI];
    //}

    //scalar mixtureCp = 0.0;
    //scalar mixtureHs = 0.0;

    forAll(Y_,specieI)
    {
        mixture   += Y[specieI]/specieThermo_[specieI].W()*specieThermo_[specieI];
        //mixtureCp += specieThermo_[specieI].Cp(p,TInit);
        //mixtureHs += specieThermo_[specieI].He(p,TInit);
    }


     TResult = TInit;
     scalar Ttol = TInit*1e-4;
     label iter = 0;



     //Newton iteration
     do
     {

         TInit = TResult;


         //TResult = TInit - (mixture.Hs(TInit) - h)/mixture.Cp(TInit);

         TResult = TInit - (mixture.HE(p,TInit) - h)/mixture.Cp(p,TInit);
         //TResult = TInit - (mixtureHs - h)/mixtureCp;



         //clip value since calling H(T) and Cp(T) out of range is a critical error
         TResult=max(TResult,_EDC_TMin_);
         TResult=min(TResult,_EDC_TMax_);

         iter++;

     } while ((mag(TResult - TInit) > Ttol) && (iter < 100));


     if (iter >= 100)
     {
        TResult = TInit;
        Info << "Max number of iterations reached; setting TResult to " << TInit << endl;
     }


    // calculate mean mole fraction in fine structures
    scalar MMean = 0.;
    forAll(Y_,specieI)
    {
        MMean += Y[specieI]/specieThermo_[specieI].W();
    }
    MMean = 1./MMean;

    rho= p*MMean/ (TResult*mixture.RR);


}

template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_derivative
(Foam::scalarField&  YStarCell, Foam::scalarField& dYdtCell) const
{


    //Info << "edcPSR/EDC_Derivative() ... "<< endl;
    //which cell are we solving for?
    scalar cellI = YStarCell[nSpecie_];

    // surrounding fluid mass fractions
    scalarField YSurrCell(nSpecie_);


    forAll(_EDC_sfYStar_,specieI)
    {
        // clip YStarCell to positive values
        YStarCell[specieI]=max(YStarCell[specieI],0.0);
        // calculate surrounding's mass fractions
        YSurrCell[specieI]=(Y_[specieI].internalField()[cellI] - YStarCell[specieI] * _EDC_gammaStar_[cellI])/ ( 1.- _EDC_gammaStar_[cellI]);
    }

    // assume hStar = hBar !!
    //scalar hStarCell=this->thermo().hs().internalField()[cellI]; // old version OF211


    // scalar hStarCell=this->thermo().he().internalField()[cellI];
    scalar hStarCell=this->thermo().he().internalField()[cellI];



    // pStar = pBar !!
    scalar pCell=this->thermo().p().internalField()[cellI];


    //update temperature
    scalar TStarCell = 0.;
    scalar rhoStarCell = 0.;
    //Info << "edcPSR/EDC_iterateTAndRho() ... "<< endl;
    EDC_iterateTAndRho(TStarCell, rhoStarCell, hStarCell, pCell , YStarCell, this->thermo().T().internalField()[cellI]);
    //Info << "edcPSR/EDC_iterateTAndRho() ...f... "<< endl;
    // now calculate fine structure concentrations with rhoStar
    scalarField CStarCell(nSpecie_);
    //initialize YStarCell and YSurrCell from y
    forAll(_EDC_sfYStar_,specieI)
    {
        CStarCell[specieI]=max(YStarCell[specieI]*rhoStarCell/specieThermo_[specieI].W(),1e-20);
    }
    //calculate chemistry kinetics

    scalarField dCdtCell(nSpecie_);

    //Info << "edcPSR/omega() ... "<< endl;
    dCdtCell = omega(CStarCell,TStarCell,pCell);
    //Info << "edcPSR/omega() f... "<< endl;

    //convert concentrations to mass fractions
    // dYdt is chemical kinetics plus m* (Y°-Y*)


    forAll(_EDC_sfYStar_,specieI)
    {
        dYdtCell[specieI]=dCdtCell[specieI]/rhoStarCell*specieThermo_[specieI].W() +  _EDC_tau_[cellI] * (YSurrCell[specieI]-YStarCell[specieI]) ;
    }



    dYdtCell[nSpecie_]= 0.; // we are not solving for the cell index

    //Info << "edcPSR/EDC_Derivative() f ... "<< endl;


}


template<class CompType, class ThermoType>
void  Foam::chemistryModel<CompType, ThermoType>::EDC_writeCombustionTerms()
{
    //output only if it is write time
    if ( this->mesh().time().outputTime())
    {

        EDC_writeScalarField(_EDC_ReTau_,"ReTau");

	/*	
        EDC_writeScalarField(_EDC_tau_,"Tau");
        EDC_writeScalarField(_EDC_edcFactor_,"edcFactor");
        EDC_writeScalarField(_EDC_gammaStar_,"gammaStar");
        EDC_writeScalarField(_EDC_mDotStar_,"edcMdotStar");

	*/


        EDC_calcMixtureFraction();


	/*
         const scalarField rho = this->thermo().rho();


         scalar TStarCell = 0.;
         scalar rhoStarCell = 0.;
         scalar mixtuteCpStarCell = 0.;
         scalar mixtuteHsStarCell = 0.;
         scalar mixtuteRRStarCell = 0.;


         scalarField YStarCell(nSpecie_);

         scalarField TStar(rho.size());
         scalarField mixtureCpStar(rho.size());
         scalarField mixtureHsStar(rho.size());
         scalarField mixtureRRStar(rho.size());


         forAll(TStar,cellI)
         {
             dimensionedScalar hCell=this->thermo().he().internalField()[cellI];
             scalar pCell=this->thermo().p().internalField()[cellI];
             scalar TCell=this->thermo().T().internalField()[cellI];

             forAll(Y_,specieI)
             {
                 YStarCell[specieI]=_EDC_sfYStar_[specieI][cellI];
             }
             EDC_iterateTAndRho(TStarCell, rhoStarCell, this->thermo().he().internalField()[cellI],  pCell , YStarCell, TCell, mixtuteCpStarCell, mixtuteRRStarCell, mixtuteHsStarCell );

             TStar[cellI]=TStarCell;

             mixtureCpStar[cellI] = mixtuteCpStarCell;
             mixtureRRStar[cellI] = mixtuteRRStarCell;
             mixtureHsStar[cellI] = mixtuteHsStarCell;

         }

         EDC_writeScalarField(TStar,"EDC_TStar");
         EDC_writeScalarField(mixtureCpStar,"EDC_mixtureCpStar");
         EDC_writeScalarField(mixtureRRStar,"EDC_mixtureRRStar");
         EDC_writeScalarField(mixtureHsStar,"EDC_mixtureHsStar");

         //const basicMultiComponentMixture& mixture = dynamic_cast<const basicMultiComponentMixture&>(this->thermo());
         // mixture.species()[specieI];


         forAll(_EDC_RR_, specieI)
         {
             // EDC_writeScalarField(RR_[specieI],"EDC_RR." + specieThermo_[specieI].name());
             // EDC_writeScalarField(_EDC_vfYStar_[specieI],"EDC_YStar." + specieThermo_[specieI].name());

              //EDC_writeScalarField(RR_[specieI],"EDC_RR." + mixture.species()[specieI]);


              EDC_writeScalarField(RR_[specieI],"RR." + this->Y_[specieI].name());
              EDC_writeScalarField(_EDC_RR_[specieI],"EDC_RR." + this->Y_[specieI].name());
              EDC_writeScalarField(_EDC_sfYStar_[specieI],"EDC_YStar." + this->Y_[specieI].name());


         }




         const scalarField h = this->thermo().he();
         EDC_writeScalarField(h,"EDC_h");

	*/


    } // end if 
}


template<class CompType, class ThermoType>
void  Foam::chemistryModel<CompType, ThermoType>::EDC_writeScalarField(scalarField X, word fieldName)
{
    //output only if it is write time
    //if ( this->mesh().time().outputTime())
    {
        // as paraView only undestands volScalarFields, make a volScalarField out of the scalarField and write it
        volScalarField volX
        (
                IOobject
                (
                    fieldName,
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
        );

        volX.internalField() = X;
        volX.correctBoundaryConditions();
        volX.write();
    }
}


template<class CompType, class ThermoType>
void  Foam::chemistryModel<CompType, ThermoType>::EDC_calcMixtureFraction()
{


    // Update 11 September 2014
    // for mixture fraction calculation
    // Global elements C,H,N,O, AR counts in a mixture



     // write mixture fraction
    const scalarField rho = this->thermo().rho();
    scalarField mixtureFraction(rho.size());

    scalarField YC(rho.size());
    scalarField YH(rho.size());

    const scalar Y1H = 0.0393;
    const scalar Y1C = 0.1170;
    const scalar Y2H = 0.0007;
    const scalar Y2C = 0.0000;

    const scalar WTC = 12.01115;
    const scalar WTH = 1.00797;


    scalar cellYH = 0.0;
    scalar cellYC = 0.0;

    //forAll(mixtureFraction,cellI)
    for (unsigned int cellI = 0; cellI< rho.size(); cellI++)
        {

         mixtureFraction[cellI] = cellYH = cellYC =  0;


         for (unsigned int k = 0; k< _EDC_specieAtomsTable_.size(); k++)
             {
        //if  (el_data[k].symbol == 'C')
        //if  (el_data[k].symbol == 'H')

        cellYC +=  _EDC_specieAtomsTable_[k]['C'] * Y_[k].internalField()[cellI] * WTC / specieThermo_[k].W();
        cellYH +=  _EDC_specieAtomsTable_[k]['H'] * Y_[k].internalField()[cellI] * WTH / specieThermo_[k].W();

             }


         YC[cellI] = cellYC ;
         YH[cellI] = cellYH ;

         const scalar denominator = 0.5 * (Y1H -Y2H) / WTH + 2.0 * (Y1C - Y2C) / WTC ;

             mixtureFraction[cellI] =  ( 0.5 * (cellYH -Y2H) / WTH + 2.0 * (cellYC - Y2C) / WTC ) / denominator  ;

         if (mixtureFraction[cellI] < 0)
            mixtureFraction[cellI] = 0;
    }



    //EDC_writeScalarField(mixtureFraction,"mixtureFraction");
    EDC_writeScalarField(YC,"YCelemental");
    EDC_writeScalarField(YH,"YHelemental");


}



template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_LES_solve()
{


    //const scalar Cmu  = 0.09;
    const scalar X    = 1.00;

    Info << "EDC/Perfectly Stirred Reactor for LES" << endl;

    // get flow parameters: rho, k, epsilon, nu
    const scalarField rho = this->thermo().rho();

    //const scalarField& k = this->db().objectRegistry::lookupObject<volScalarField>("kSgs").internalField();
    //const scalarField& epsilon = this->db().objectRegistry::lookupObject<volScalarField>("eSgs").internalField();
    const scalarField& tauStar = this->db().objectRegistry::lookupObject<volScalarField>("tauStar").internalField();

    //const scalarField nu = this->thermo().mu().internalField()/rho;

    // calculate EDC model variables: gammaStar and mDotStar
    const scalarField gammaL = _EDC_CD1_ * this->db().objectRegistry::lookupObject<volScalarField>("gammaL").internalField();

    _EDC_gammaStar_= pow(gammaL,2.);
    _EDC_gammaStar_= min (_EDC_gammaStarClipFactor_, _EDC_gammaStar_);

    _EDC_mDotStar_ = _EDC_CD2_ * pow (tauStar, -1.0) ;

    _EDC_tau_  =  _EDC_CD3_ * _EDC_mDotStar_;

    // calculate reaction rate coefficient
    _EDC_edcFactor_ = rho * _EDC_gammaStar_ * _EDC_mDotStar_ * X  / (1.0 - _EDC_gammaStar_*X);
    //edcFactor_ = rho * gammaStar_ * mDotStar_ * X ;

    // calculate fine structure composition according to fast chemistry approach
    Info << "Integrating PSR equations from 0 to steady state ... ";
    EDC_updateYStar();
    Info << "done" << endl;
    //transform reaction rate to un-normalized specie space



    forAll(rho, celli)
    {

        for (label i=0; i<nSpecie_; i++)
        {
            //RR_[i][celli] = (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_vfYStar_[i][celli]);
                  RR_[i][celli] =  (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_sfYStar_[i][celli]);
             _EDC_RR_[i][celli] =  (-1.0) * _EDC_edcFactor_[celli] * (Y_[i][celli] - _EDC_sfYStar_[i][celli]);


              //Info << "specie\t" << i << "\tRR:" << RR_[i][celli]*1e+10 << "\tEDCRR:" << _EDC_RR_[i][celli]*1e+10  <<   endl;

        }
    }




}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::EDC_iterateTAndRho
(Foam::scalar& TResult, Foam::scalar& rho, Foam::scalar h, Foam::scalar p, Foam::scalarField& Y, Foam::scalar TInit,
    Foam::scalar& mixtureCp, Foam::scalar& mixtureRR, Foam::scalar& mixtureHs) const
{

     ThermoType mixture(0.0*specieThermo_[0]);
     // use size from global variable, so that YStarCell (which contains the cell index as Y_.size+1)

     //scalar mixtureCp;
     //scalar

     forAll(Y_,specieI)
     {
         mixture   += Y[specieI]/specieThermo_[specieI].W()*specieThermo_[specieI];
         //mixtureCp += specieThermo_[specieI].Cp(p,TInit);
         //mixtureHs += specieThermo_[specieI].Hs(p,TInit);
     }




     TResult = TInit;
     scalar Ttol = TInit*1e-4;
     label iter = 0;

     //Newton iteration
     do
     {

         TInit = TResult;
         TResult = TInit - (mixture.HE(p,TInit) - h)/mixture.Cp(p,TInit);
         //TResult = TInit - (mixtureHs - h)/mixtureCp;

         //clip value since calling H(T) and Cp(T) out of range is a critical error
         TResult=max(TResult,_EDC_TMin_);
         TResult=min(TResult,_EDC_TMax_);

         iter++;

     } while ((mag(TResult - TInit) > Ttol) && (iter < 100));


     if (iter >= 100)
     {
        TResult = TInit;
        Info << "Max number of iterations reached; setting TResult to " << TInit << endl;
     }


    // calculate mean mole fraction in fine structures
    scalar MMean=0.;
    forAll(Y_,specieI)
    {
        MMean += Y[specieI]/specieThermo_[specieI].W();
    }
    MMean = 1./MMean;

    rho= p*MMean/ (TResult*mixture.RR);


    /*
    forAll(Y_,specieI)
    {
        //mixture   += Y[specieI]/specieThermo_[specieI].W()*specieThermo_[specieI];
        mixtureCp += specieThermo_[specieI].Cp(p,TResult);
        mixtureHs += specieThermo_[specieI].Hs(p,TResult);
    }*/

    mixtureCp = mixture.Cp(p,TResult);
    mixtureHs = mixture.HE(p,TResult);
    mixtureRR = mixture.RR;


}

// ************************************************************************* //

template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    const scalar t0,
    const scalar deltaT
)
{
    Info << endl << endl << endl << "************** SOLVE FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
    return 0;

}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistryModel<CompType, ThermoType>::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{


     Info << endl << endl << endl << "************** SOLVE FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
     return 0;

}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::updateConcsInReactionI
(
    const label index,
    const scalar dt,
    const scalar omeg,
    const scalar p,
    const scalar T,
    scalarField& c
) const
{
     Info << endl << endl << endl << "************** SOLVE FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
}


template<class CompType, class ThermoType>
void Foam::chemistryModel<CompType, ThermoType>::updateRRInReactionI
(
    const label index,
    const scalar pr,
    const scalar pf,
    const scalar corr,
    const label lRef,
    const label rRef,
    const scalar p,
    const scalar T,
    simpleMatrix<scalar>& RR
) const
{
     Info << endl << endl << endl << "************** SOLVE FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
}






