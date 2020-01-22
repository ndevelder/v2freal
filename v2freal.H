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

Class
    Foam::incompressible::RASModels::v2freal

Description
    Implementation of the v2freal turbulence model for incompressible
    flows.

    Turbulence model described in:
    @verbatim
        Perot, J. Blair
        "Turbulence modeling using body force potentials"
        Physic of Fluids,
        1999
    @endverbatim

    The default model coefficients correspond to the following:
    @verbatim
        v2frealCoeffs
        {

        }
    @endverbatim

SourceFiles
    v2freal.C

\*---------------------------------------------------------------------------*/

#ifndef v2freal_H
#define v2freal_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class v2freal Declaration
\*---------------------------------------------------------------------------*/

class v2freal
:
    public RASModel
{
    // Private data

        // Model coefficients
            dimensionedScalar cEp1_;
            dimensionedScalar cEp2con_;
            dimensionedScalar cP1_;
            dimensionedScalar cP2_;
            dimensionedScalar cP3_;
            dimensionedScalar cP4_;
			dimensionedScalar cL1_;
			dimensionedScalar cL2_;
			dimensionedScalar cLim_;
            dimensionedScalar cMu_;
			dimensionedScalar betaK_;
			dimensionedScalar nutRatMax_;
            dimensionedScalar sigmaKInit_;
            dimensionedScalar sigmaEpsInit_;
            dimensionedScalar sigmaV2Init_;
            word              solveK_;
            word              solveEps_;
            word              solveV2_;
            word              solveNut_;
            word              eqnSigmaK_;
            word              eqnSigmaEps_;
            word              eqnSigmaV2_;
			word              realNut_;
			word              durbinReal_;
			word              debugWrite_;

        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
            wallDist y_;

        // Fields

            volScalarField k_;
            volVectorField gradk_;
            volScalarField epsilon_;
            volScalarField nut_;
            volScalarField nutNorm_;
            volScalarField v2_;
            volVectorField vorticity_;
            volTensorField uGrad_;
			volScalarField magS_;
			volScalarField kSqrt_;
			volVectorField gradkSqrt_;
			volScalarField sigmaK_;
			volScalarField sigmaEps_;
			volScalarField sigmaV2_;
			volScalarField cEp2_;
			volScalarField f_;
			volScalarField G_;







    // Private member functions

    //- Return time scale, Ts using epsilon
    tmp<volScalarField> Ts() const;
    tmp<volScalarField> Ls() const;
	tmp<volScalarField> Tsdurbin() const;
	tmp<volScalarField> Lsdurbin() const;
	
	tmp<volScalarField> davidsonCorrectNut
        (
            const tmp<volScalarField>& value
		) const;
    
    //- Return time scale, TsEh using epsHat
    //tmp<volScalarField> TsEh() const;


public:

    //- Runtime type information
    TypeName("v2freal");


    // Constructors

        //- Construct from components
        v2freal
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor

        virtual ~v2freal()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }


        //- Return the local turbulent reynolds number
        tmp<volScalarField> reTau() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("reTau", (k_*k_ + k0_*k0_)/(nu()*(epsilon_ + epsilonSmall_)))
            );
        }


        //- Return cEp2
        tmp<volScalarField> cEp2() const
        {
                return cEp2_;
        }

		//- Return the effective diffusivity for epsilon
        tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("nuEff", nut() + nu())
            );
        }
		
        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_*sigmaK_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_*sigmaEps_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> Dv2Eff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("Dv2Eff", nut_*sigmaV2_ + nu())
            );
        }


        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DEff", nut_ + nu())
            );
        }

        //- Return the Reynolds number
        tmp<volScalarField> tpReynolds() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("tpReynolds", k_*k_/(epsilon_*nu()))
            );
        }

        //- Return the vorticity
        virtual tmp<volVectorField> vorticity() const
        {
            return vorticity_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
		
        //- Return f
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }
		
        //- Return v2
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

		//- Return the turbulence kinetic energy
        virtual tmp<volVectorField> U() const
        {
            return U_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
