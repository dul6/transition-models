/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::RASModels::kOmegaSSTGam

Group
    grpLESTurbulence

Description
    A One-Equation Local Correlation-Based Transition Model
    based on the k-omega-SST RAS model.

    References:
    \verbatim
	Menter, F.R., Smirnov, P.E., Liu, T. et al. (2015)
	A One-Equation Local Correlation-Based Transition Model.
	Flow Turbulence Combust 95, 583–619.
    \endverbatim

SourceFiles
    kOmegaSSTGam.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTGam_H
#define kOmegaSSTGam_H

#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class kOmegaSSTGam Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTGam
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        kOmegaSSTGam(const kOmegaSSTGam&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTGam&) = delete;


protected:

    // Protected data

        // Model constants

            dimensionedScalar ca2_;

            dimensionedScalar ce2_;

        // Fields

            //- Intermittency
            volScalarField gammaInt_;

            //- Effective intermittency
            volScalarField::Internal gammaIntEff_;

            const volVectorField& n_;

            volScalarField::Internal P_klim_;

    // Protected Member Functions

        //- Modified form of the k-omega SST F1 function
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Modified form of the k-omega SST epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Empirical correlation for critical Reynolds number where the
        //  intermittency first starts to increase in the boundary layer
        tmp<volScalarField::Internal> ReThetac
        (
            const volScalarField::Internal& dVdy,
            const volScalarField::Internal& nu
        ) const;

        //- Transition onset location control function
        tmp<volScalarField::Internal> Fonset
        (
            const volScalarField::Internal& Rev,
            const volScalarField::Internal& ReThetac,
            const volScalarField::Internal& RT
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correctReThetatGammaInt();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTGam");


    // Constructors

        //- Construct from components
        kOmegaSSTGam
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTGam() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to intermittency
        const volScalarField& gammaInt() const
        {
            return gammaInt_;
        }

        //- Return the effective diffusivity for intermittency
        tmp<volScalarField> DgammaIntEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DgammaIntEff",
                    this->nut_ + this->nu()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTGam.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
