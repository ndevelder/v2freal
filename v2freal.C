/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "v2freal.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hello

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(v2freal, 0);
addToRunTimeSelectionTable(RASModel, v2freal, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> v2freal::Ts() const
{ 
    return max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)));
}


tmp<volScalarField> v2freal::Ls() const
{	
	return cL1_*max(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
}

tmp<volScalarField> v2freal::Tsdurbin() const
{ 
    return min((max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)))), (cLim_/sqrt(3.0))*k_/((v2_+k0_)*cMu_*magS_));
}


tmp<volScalarField> v2freal::Lsdurbin() const
{	
	return cL1_*max(min(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_), (1.0/sqrt(3.0))*pow(k_+k0_, 1.5)/((v2_+k0_)*cMu_*magS_)),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2freal::v2freal
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.4
        )
    ),
    cEp2con_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.9
        )
    ),
    cP1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            1.4
        )
    ),
    cP2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.3
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cL1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL1",
            coeffDict_,
            0.23
        )
    ),
    cL2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL2",
            coeffDict_,
            70.0
        )
    ),
    cLim_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cLim",
            coeffDict_,
            1.0
        )
    ),
    cMu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.22
        )
    ),
    betaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaK",
            coeffDict_,
            0.09
        )
    ),
	nutRatMax_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutRatMax",
            coeffDict_,
            1.0e5
        )
    ),
    sigmaKInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),
    sigmaEpsInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.769
        )
    ),
    sigmaV2Init_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaV2Init",
            coeffDict_,
            0.83
        )
    ),

   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solveV2_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaV2_
   (
       coeffDict_.lookup("eqnSigmaV2")
   ),
   
   realNut_
   (
       coeffDict_.lookup("realNut")
   ),
   
   durbinReal_
   (
       coeffDict_.lookup("durbinReal")
   ),

   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   y_
   (
   mesh_
   ),
    
	k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(k_))
    ),
    
	epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_/max(nut_))
    ),
    
	v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
    
	vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::curl(U_))
    ),
    
	uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::grad(U_))
    ),
    
	magS_
    (
        IOobject
        (
            "magS",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(fvc::grad(U_))
    ),
	
	kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(k_))
    ),
    
	gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(kSqrt_))
    ),
    
	sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaKInit_)
    ),
    
	sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaEpsInit_)
    ),
    
	sigmaV2_
    (
        IOobject
        (
            "sigmaV2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaV2Init_)
    ),
    
	cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
	
	f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
	
	G_
    (
        IOobject
        (
            "RASModel::G",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nut_*2*magSqr(dev(symm(uGrad_)))
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		nut_ = cMu_*v2_*Ts();	
        nut_ = min(nut_,nutRatMax_*nu());        
        nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), ROOTVSMALL));       
    }
	
    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solveV2 is: " <<solveV2_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> v2freal::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> v2freal::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> v2freal::divDevReff() const
{
    return
    (
      - fvm::laplacian(nuEff(), U_)
      - fvc::div(nuEff()*dev(fvc::grad(U_)().T()))
	);
}


bool v2freal::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		betaK_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaV2Init_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void v2freal::correct()
{

    //**********************************************//	
    // Bounding values not already defined by model
    //**********************************************//
	
	const dimensionedScalar nut0("minNut", nut_.dimensions(), ROOTVSMALL);
	const dimensionedScalar v20("minv2", v2_.dimensions(), ROOTVSMALL);
	const dimensionedScalar f0("fMin", f_.dimensions(), 0.0);
	const dimensionedScalar L0("lMin", dimensionSet(0,1,0,0,0,0,0), ROOTVSMALL);

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, k0_);
        bound(epsilon_, epsilonSmall_);
		bound(v2_,v20);
		bound(nut_,nut0);
    }
	
	
    RASModel::correct();

	
    if (!turbulence_)
    {
        return;
    }
	
	
	word sMMdebug = runTime_.controlDict().lookup("showMaxMin");


    //*************************************//	
    // Vorticity and Gradient
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);
	



    //*************************************//	
    // K Production
    //*************************************//
 
	const volScalarField S2 = 2*magSqr(dev(symm(uGrad_)));
	magS_ = sqrt(S2);
	volScalarField G("RASModel::G", nut_*S2); 
	G_ = nut_*S2;
	volScalarField GdK("GdK", G/(k_ + k0_));
	
	
    //*************************************//	
    // Length and Time Scales
    //*************************************//	
    
	volScalarField L("Length",Ls());
	volScalarField T("Time",Ts());
	
	volScalarField L_real("Length_Real",Lsdurbin());
	volScalarField T_real("Time_Real",Tsdurbin());
	
    const dimensionedScalar N("N", dimless, 6.0);	
	const volScalarField L2("Lsqr",sqr(L));
	const volScalarField L2_real("Lsqr_real",sqr(L_real));
	const volScalarField v2k("v2OverK",v2_/(k_ + k0_));	
			
	Info << "Made it past L2" << endl;
	
	
	// Make BCs function
	volScalarField tpProd("tpProd", GdK);
	volScalarField tpphi("tpphi", v2k);
	volScalarField tpphiSqrt("tpphiSqrt", sqrt(v2k));
	volVectorField tppsi("tppsi", nut_*vorticity_);

	const volScalarField pOD = G/epsilon_; 
	
	

    //*************************************//
    //Dissipation equation
    //*************************************//	
	
	volScalarField cEp1eqn("cEp1eqn",1.4*(1.0 + 0.045*min(sqrt(k_/v2_), scalar(100.0))));


    tmp<fvScalarMatrix> epsEqn  
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        cEp1eqn*G/T_real
		- fvm::Sp(cEp2con_/T_real, epsilon_)
    );

    if(solveEps_ == "true")
    {
    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,epsilonSmall_);
    }
	
	
	
	
    //*************************************//
    // Turbulent kinetic energy equation
    //*************************************//
    
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/k_,k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
    }
	
	

    
	
    //*************************************//
	// Update K-related fields
    //*************************************//
    	
    kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(ROOTVSMALL)));

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
    

	// Pressure strain terms	
	const volScalarField slowPS
    (
        "v2freal::slowPS",
        (cP1_ - 1.0)*((2.0/3.0) - v2_/(k_+k0_))/T	
	);
	
	const volScalarField fastPS
    (
        "v2freal::fastPS",
		cP2_*GdK
	);
	
	const volScalarField fwall
    (
        "v2freal::fwall",
		5.0*(v2_/(k_+k0_))/T	
	);
	
	
	
	// Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2_real, f_)
      + slowPS/L2_real 
	  + fastPS/L2_real
	  + fwall/L2_real
    );
 
    fEqn().relax();
    solve(fEqn);



    tmp<fvScalarMatrix> v2Eqn
    (	  
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
	  + fvm::SuSp(-fvc::div(phi_), v2_)
      - fvm::laplacian(Dv2Eff(), v2_)
      ==
        min(k_*f_, k_*(slowPS + fastPS + fwall))
	  - fvm::Sp(6.0*epsilon_/k_, v2_)
    );
	

    v2Eqn().relax();
    solve(v2Eqn);
	bound(v2_,v20);
	
	
	
	
    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
 
        nut_ = cMu_*v2_*T_real;
		
		// if(realNut_ == "true")
		// {
			// nut_ = min(nut_,cLim_*k_/(sqrt(6.0)*magS_));
			// Info << "Using realizability constraint" << endl;
		// }
		
		nut_ = min(nut_,nutRatMax_*nu());
        nut_.correctBoundaryConditions();		
    }
	
	
	
    //*************************************//   
    // Output some max values
    //*************************************//
	
	if(sMMdebug == "true")
	{
    
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	
	Info<< "Max T: " << gMax(T) << " Min T: " << gMin(T) << endl;
	Info<< "Max L: " << gMax(L) << " Min L: " << gMin(L) << endl;
	Info << "Max cEp1: " << max(cEp1eqn) << " Min cEp1: " << min(cEp1eqn) << endl; 
	Info<< "Max f: " << gMax(f_) << " Min f: " << gMin(f_) << endl;
    Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Epsilon: " << gMax(epsilon_) <<endl;
    Info<< "Max v2: " << gMax(v2_) << " Min v2: " << gMin(v2_) << " Max G: " << gMax(G) <<endl;
    Info<< "Max uTauSquared: " << gMax(uTauSquared) << " Max vorticity: " << gMax(vorticity_) << endl;
	
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
