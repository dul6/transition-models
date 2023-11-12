/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "kOmegaSSTGamDES.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTGamDES<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    const volScalarField Ry(this->y_*sqrt(this->k_)/this->nu());
    const volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(kOmegaSSTDES<BasicTurbulenceModel>::F1(CDkOmega), F3);
}
 

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTGamDES<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return gammaIntEff_*kOmegaSSTDES<BasicTurbulenceModel>::Pk(G) + P_klim_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTGamDES<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return
        max(gammaIntEff_, scalar(0.1))
       *kOmegaSSTDES<BasicTurbulenceModel>::epsilonByk(F1, gradU);
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
kOmegaSSTGamDES<BasicTurbulenceModel>::ReThetac
(
        const volScalarField::Internal& dVdy,
        const volScalarField::Internal& nu
) const
{
    tmp<volScalarField::Internal> tReThetac
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::groupName("ReThetac", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimless
        )
    );
    volScalarField::Internal& ReThetac = tReThetac.ref();

    const volScalarField& k = this->k_;
    const volScalarField::Internal& omega = this->omega_();
    const volScalarField::Internal& y = this->y_();

    forAll(ReThetac, celli)
    {
        const scalar TuL
        (
            min(100*sqrt((2.0/3.0)*k[celli])/(omega[celli]*y[celli]), scalar(100))
        );

        scalar lambdaL
        (
            -0.00757*sqr(y[celli])*dVdy[celli]/nu[celli] + 0.0128
        );

        lambdaL = min(max(lambdaL, -1.0), 1.0);

        scalar Fpg =
             lambdaL >= 0
            ?
             min(1 + 14.68*lambdaL, 1.5)
            :
             min(1 - 7.34*lambdaL, 3);
        Fpg = max(Fpg, 0);

        ReThetac[celli] = 100 + 1000*exp(-TuL*Fpg);
    }

    return tReThetac;
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTGamDES<BasicTurbulenceModel>::Fonset
(
    const volScalarField::Internal& Rev,
    const volScalarField::Internal& ReThetac,
    const volScalarField::Internal& RT
) const
{
    const volScalarField::Internal Fonset1(Rev/(2.2*ReThetac));

    const volScalarField::Internal Fonset2
    (
        min(Fonset1, scalar(2))
    );

    const volScalarField::Internal Fonset3(max(1 - pow3(RT/3.5), scalar(0)));

    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal
        (
            IOobject::groupName("Fonset", this->alphaRhoPhi_.group()),
            max(Fonset2 - Fonset3, scalar(0))
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTGamDES<BasicTurbulenceModel>::kOmegaSSTGamDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    ca2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50
        )
    ),

    gammaInt_
    (
        IOobject
        (
            IOobject::groupName("gammaInt", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    gammaIntEff_
    (
        IOobject
        (
            IOobject::groupName("gammaIntEff", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    n_(wallDist::New(this->mesh_).n()),

    P_klim_
    (
        IOobject
        (
            IOobject::groupName("P_klim", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0, 0, 0), Zero)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTGamDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        ca2_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmegaSSTGamDES<BasicTurbulenceModel>::correctReThetatGammaInt()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y_();
    fv::options& fvOptions(fv::options::New(this->mesh_));
    const volScalarField& nut = this->nut_;

    // Fields derived from the velocity gradient
    const volScalarField::Internal dVdy(fvc::grad(U & n_) & n_);
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    tgradU.clear();

    const volScalarField::Internal ReThetac(this->ReThetac(dVdy, nu));
    const volScalarField::Internal Rev(sqr(y)*S/nu);
    const volScalarField::Internal RT(k()/(nu*omega()));
    const volScalarField::Internal F_limon(min(max(Rev/(2.2 * 1100.0) -1.0, 0.0), 3.0));


    {
        const volScalarField::Internal Pgamma
        (
            alpha()*rho()
           *100*S*gammaInt_()*Fonset(Rev, ReThetac, RT)
        );

        const volScalarField::Internal Fturb(exp(-pow4(0.5*RT)));

        const volScalarField::Internal Egamma
        (
            alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
        );

        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)
        ==
            Pgamma - fvm::Sp(Pgamma, gammaInt_)
          + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
          + fvOptions(alpha, rho, gammaInt_)
        );

        gammaIntEqn.ref().relax();
        fvOptions.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvOptions.correct(gammaInt_);
        bound(gammaInt_, 0);
    }

    dimensionedScalar zMin("zMin", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0);
    const volScalarField::Internal pkk(5*max(gammaInt_ - 0.2, scalar(0.0))*(1 - gammaInt_));
    P_klim_= F_limon*pkk*Omega*S*max(3*nu - nut, zMin);//(5*max(gammaInt_ - 0.2, scalar(0.0))*(1 - gammaInt_)*F_limon*max(3*nu - nut, zMin)*Omega*S);
    gammaIntEff_ = gammaInt_();
}


template<class BasicTurbulenceModel>
void kOmegaSSTGamDES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Correct k and omega
    kOmegaSSTDES<BasicTurbulenceModel>::correct();

    // Correct ReThetat and gammaInt
    correctReThetatGammaInt();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
