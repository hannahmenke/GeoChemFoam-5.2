/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

\*---------------------------------------------------------------------------*/

#include "basicTwoPhaseMultiComponentTransportMixture.H"
#include "upwind.H"
#include "downwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicTwoPhaseMultiComponentTransportMixture, 0);
    defineRunTimeSelectionTable(basicTwoPhaseMultiComponentTransportMixture, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicTwoPhaseMultiComponentTransportMixture::basicTwoPhaseMultiComponentTransportMixture
(
    const fvMesh& mesh,
    const volScalarField& alpha1
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoPhysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    basicMultiComponentMixture(*this, this->subDict("solutionSpecies").toc(), mesh),
    mesh_(mesh),
    alpha1_(alpha1),
    cYi_
    (
        alpha1.mesh().solverDict("Yi").get<scalar>("cYi")
    ),
    D1Y_(this->subDict("solutionSpecies").toc().size()),
    D2Y_(this->subDict("solutionSpecies").toc().size()),
    HY_(this->subDict("solutionSpecies").toc().size()),
    Mw_(this->subDict("solutionSpecies").toc().size()),
    grotthussActive_(this->subDict("solutionSpecies").toc().size(), false),
    grotthussDvehicle_(this->subDict("solutionSpecies").toc().size(), 0.0),
    grotthussDstructural_(this->subDict("solutionSpecies").toc().size(), 0.0),
    grotthussSaturationSuppression_(false),
    grotthussSaturationExponent_(1.0),
    grotthussMinWaterSaturation_(1e-3),
    grotthussPHScaling_(false),
    grotthussHplusSpecies_("H+"),
    grotthussPHReference_(7.0),
    grotthussPHWidth_(1.0),
    grotthussPHMinFactor_(0.0),
    grotthussIonicStrengthScaling_(false),
    grotthussIonicStrengthSpecies_(),
    grotthussIonicStrengthCharges_(),
    grotthussIonicStrengthReference_(0.01),
    grotthussIonicStrengthExponent_(0.0),
    grotthussIonicStrengthMin_(1e-12),
    grotthussIonicStrengthMinFactor_(0.0),
    phiD_
	(
		IOobject
		(
			"phiD",
            alpha1_.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("phiD", dimMass/dimTime, 0.0)
	),
    Mflux_
	(
		IOobject
		(
			"Mflux",
            alpha1_.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("MFlux", dimMass/dimVolume/dimTime, 0.0)
	)
{

    //******************************************************************//
    // scheme for phiH
    //******************************************************************//

	OStringStream buf;

	const tokenList& tokens = mesh_.divScheme("div(phiH,Yi)");

	label len = tokens.size();

	forAll (tokens, tok)
	{
        buf << mesh_.divScheme("div(phiH,Yi)")[tok];

		if (--len)
		{
	        buf << ' ';
	    }
	}
  
    phiHScheme_ = buf.str();

	forAll(species_, i)
	{
	    const dictionary& solutionSpeciesDict = this->subDict("solutionSpecies");
		const dictionary& subdict = solutionSpeciesDict.subDict(species_[i]);
		D1Y_.set
		(
			i,
			new dimensionedScalar("D1",subdict)
		);
		D2Y_.set
		(
			i,
			new dimensionedScalar("D2",subdict)
		);
		HY_.set
		(
			i,
			new dimensionedScalar("H",subdict)
		);
		Mw_.set
		(
			i,
			new dimensionedScalar("Mw",subdict)
		);
	}

    if (this->found("grotthussTransport"))
    {
        const dictionary& grotthussDict = this->subDict("grotthussTransport");
        const Switch enabled(grotthussDict.lookupOrDefault<Switch>("enabled", false));

        if (enabled)
        {
            const scalar accelerationFactor =
                grotthussDict.lookupOrDefault<scalar>("accelerationFactor", 1.0);
            grotthussSaturationSuppression_ =
                grotthussDict.lookupOrDefault<Switch>("useWaterSaturationSuppression", false);
            grotthussSaturationExponent_ =
                grotthussDict.lookupOrDefault<scalar>("saturationExponent", 1.0);
            grotthussMinWaterSaturation_ =
                grotthussDict.lookupOrDefault<scalar>("minWaterSaturation", 1e-3);
            grotthussPHScaling_ =
                grotthussDict.lookupOrDefault<Switch>("usePHScaling", false);
            grotthussHplusSpecies_ =
                grotthussDict.lookupOrDefault<word>("hplusSpecies", "H+");
            grotthussPHReference_ =
                grotthussDict.lookupOrDefault<scalar>("pHReference", 7.0);
            grotthussPHWidth_ =
                grotthussDict.lookupOrDefault<scalar>("pHWidth", 1.0);
            grotthussPHMinFactor_ =
                grotthussDict.lookupOrDefault<scalar>("pHMinFactor", 0.0);
            grotthussIonicStrengthScaling_ =
                grotthussDict.lookupOrDefault<Switch>("useIonicStrengthScaling", false);
            grotthussIonicStrengthReference_ =
                grotthussDict.lookupOrDefault<scalar>("ionicStrengthReference", 0.01);
            grotthussIonicStrengthExponent_ =
                grotthussDict.lookupOrDefault<scalar>("ionicStrengthExponent", 0.0);
            grotthussIonicStrengthMin_ =
                grotthussDict.lookupOrDefault<scalar>("ionicStrengthMin", 1e-12);
            grotthussIonicStrengthMinFactor_ =
                grotthussDict.lookupOrDefault<scalar>("ionicStrengthMinFactor", 0.0);

            if (grotthussDict.found("ionicStrengthSpecies"))
            {
                grotthussDict.lookup("ionicStrengthSpecies") >> grotthussIonicStrengthSpecies_;
            }
            if (grotthussDict.found("ionicStrengthCharges"))
            {
                grotthussDict.lookup("ionicStrengthCharges") >> grotthussIonicStrengthCharges_;
            }

            if (grotthussMinWaterSaturation_ < 0.0 || grotthussMinWaterSaturation_ >= 1.0)
            {
                FatalErrorInFunction
                    << "grotthussTransport minWaterSaturation must be in [0, 1), got "
                    << grotthussMinWaterSaturation_ << exit(FatalError);
            }

            if (grotthussSaturationExponent_ <= 0.0)
            {
                FatalErrorInFunction
                    << "grotthussTransport saturationExponent must be positive, got "
                    << grotthussSaturationExponent_ << exit(FatalError);
            }

            if (grotthussPHScaling_ && grotthussPHWidth_ <= 0.0)
            {
                FatalErrorInFunction
                    << "grotthussTransport pHWidth must be positive, got "
                    << grotthussPHWidth_ << exit(FatalError);
            }

            if (grotthussPHMinFactor_ < 0.0 || grotthussPHMinFactor_ > 1.0)
            {
                FatalErrorInFunction
                    << "grotthussTransport pHMinFactor must be in [0, 1], got "
                    << grotthussPHMinFactor_ << exit(FatalError);
            }

            if (grotthussIonicStrengthScaling_)
            {
                if
                (
                    grotthussIonicStrengthSpecies_.size()
                 != grotthussIonicStrengthCharges_.size()
                 || grotthussIonicStrengthSpecies_.empty()
                )
                {
                    FatalErrorInFunction
                        << "grotthussTransport ionic strength scaling requires "
                        << "nonempty ionicStrengthSpecies and ionicStrengthCharges "
                        << "lists of the same length." << exit(FatalError);
                }

                if (grotthussIonicStrengthReference_ <= 0.0)
                {
                    FatalErrorInFunction
                        << "grotthussTransport ionicStrengthReference must be positive, got "
                        << grotthussIonicStrengthReference_ << exit(FatalError);
                }

                if (grotthussIonicStrengthMin_ <= 0.0)
                {
                    FatalErrorInFunction
                        << "grotthussTransport ionicStrengthMin must be positive, got "
                        << grotthussIonicStrengthMin_ << exit(FatalError);
                }

                if
                (
                    grotthussIonicStrengthMinFactor_ < 0.0
                 || grotthussIonicStrengthMinFactor_ > 1.0
                )
                {
                    FatalErrorInFunction
                        << "grotthussTransport ionicStrengthMinFactor must be in [0, 1], got "
                        << grotthussIonicStrengthMinFactor_ << exit(FatalError);
                }
            }

            wordList grotthussKeys(grotthussDict.toc());

            forAll(species_, i)
            {
                forAll(grotthussKeys, keyi)
                {
                    if (!grotthussDict.isDict(grotthussKeys[keyi]))
                    {
                        continue;
                    }

                    const dictionary& speciesDict =
                        grotthussDict.subDict(grotthussKeys[keyi]);

                    const word transportName =
                        speciesDict.lookupOrDefault<word>("name", grotthussKeys[keyi]);

                    if (transportName != species_[i])
                    {
                        continue;
                    }

                    dimensionedScalar Dvehicle("Dvehicle", speciesDict);
                    dimensionedScalar Dstructural("Dstructural", speciesDict);
                    const scalar DvehicleValue = accelerationFactor*Dvehicle.value();
                    const scalar DstructuralValue = accelerationFactor*Dstructural.value();

                    D1Y_.set
                    (
                        i,
                        new dimensionedScalar
                        (
                            "D1",
                            Dvehicle.dimensions(),
                            DvehicleValue + DstructuralValue
                        )
                    );
                    grotthussActive_[i] = true;
                    grotthussDvehicle_[i] = DvehicleValue;
                    grotthussDstructural_[i] = DstructuralValue;

                    Info<< "Grotthuss effective D1 for " << species_[i]
                        << " = " << D1Y_[i].value() << " m2/s" << nl;
                    if (grotthussSaturationSuppression_)
                    {
                        Info<< "Grotthuss structural diffusion for " << species_[i]
                            << " is scaled by water saturation with exponent "
                            << grotthussSaturationExponent_ << " and minimum saturation "
                            << grotthussMinWaterSaturation_ << nl;
                    }
                    if (grotthussPHScaling_)
                    {
                        Info<< "Grotthuss structural diffusion for " << species_[i]
                            << " is scaled by pH distance from "
                            << grotthussPHReference_ << " with pHWidth "
                            << grotthussPHWidth_ << " and pHMinFactor "
                            << grotthussPHMinFactor_ << nl;
                    }
                    if (grotthussIonicStrengthScaling_)
                    {
                        Info<< "Grotthuss structural diffusion for " << species_[i]
                            << " is scaled by ionic strength using species "
                            << grotthussIonicStrengthSpecies_
                            << " and reference "
                            << grotthussIonicStrengthReference_ << nl;
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicTwoPhaseMultiComponentTransportMixture::~basicTwoPhaseMultiComponentTransportMixture()
{}

// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //
//calculate and return the henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::phiH(const label i) const
{
    return DmY(i)*(1-HY_[i])/fvc::interpolate((alpha1_+HY_[i]*(1-alpha1_)))*fvc::snGrad(alpha1_)*mesh_.magSf();
}

//calculate and return the upwind henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::phiHUp(const label i) const
{
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();
	surfaceScalarField alphaUp = upwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    return DmY(i)/(alphaUp+(1-alphaUp)*HY_[i])*fvc::snGrad(alpha1_)*mesh_.magSf();
}

//calculate and return the downwind henry transfer flux
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture:: phiHDown(const label i) const
{
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();
	surfaceScalarField alphaDown = downwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
    return -HY_[i]*DmY(i)/(alphaDown+(1-alphaDown)*HY_[i])*fvc::snGrad(alpha1_)*mesh_.magSf();

}

//calculate compression coeeficient for species i
surfaceScalarField Foam::basicTwoPhaseMultiComponentTransportMixture::compressionCoeff(const label i)
{
    //Direction of interfacial flux
    surfaceScalarField fluxDir = fvc::snGrad(alpha1_)*mesh_.magSf();

    //upwind and downwind alpha1
	surfaceScalarField alphaUp = upwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);
	surfaceScalarField alphaDown = downwind<scalar>(mesh_,fluxDir).interpolate(alpha1_);

    //upwind and downwnd Yi
    surfaceScalarField YiUp   = upwind<scalar>(mesh_,fluxDir).interpolate(Y_[i]);
    surfaceScalarField YiDown = downwind<scalar>(mesh_,fluxDir).interpolate(Y_[i]);
    
    dimensionedScalar sgn = sign(max(alphaDown*YiDown)-max((1-alphaUp)*YiUp));

    //normal compression coefficient
    surfaceScalarField deltaYi1 = max(-max(Y_[i]),min(max(Y_[i]),(YiDown-YiUp)/(alphaDown-alphaUp+1e-4)));
    
    //standard compression coefficient
    surfaceScalarField deltaYi2 = max(-max(Y_[i]),min(max(Y_[i]),YiDown/(alphaDown+(1-alphaDown)*HY_[i])-HY_[i]*YiUp/(alphaUp+(1-alphaUp)*HY_[i])));

    return sgn*max(mag(deltaYi1),mag(deltaYi2));

}
