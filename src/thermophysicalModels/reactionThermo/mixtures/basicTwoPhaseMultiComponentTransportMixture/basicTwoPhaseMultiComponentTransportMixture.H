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

Class
    Foam::basicTwoPhaseMultiComponentTransportMixture

Description
    Foam::basicTwoPhaseMultiComponentTransportMixture

SourceFiles
    basicTwoPhaseMultiComponentTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicTwoPhaseMultiComponentTransportMixture_H
#define basicTwoPhaseMultiComponentTransportMixture_H

#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicTwoPhaseMultiComponentTransportMixture Declaration
\*---------------------------------------------------------------------------*/

class basicTwoPhaseMultiComponentTransportMixture
:
    public IOdictionary,
    public basicMultiComponentMixture
{

protected:

    // Protected data
        const fvMesh& mesh_;
        const volScalarField& alpha1_;

        scalar cYi_;

        string phiHScheme_;

        //- Diffusion and Henry's constant
        PtrList<dimensionedScalar> D1Y_;
        PtrList<dimensionedScalar> D2Y_;
        PtrList<dimensionedScalar> HY_;
        PtrList<dimensionedScalar> Mw_;
        List<bool> grotthussActive_;
        List<scalar> grotthussDvehicle_;
        List<scalar> grotthussDstructural_;
        bool grotthussSaturationSuppression_;
        scalar grotthussSaturationExponent_;
        scalar grotthussMinWaterSaturation_;
        bool grotthussPHScaling_;
        word grotthussHplusSpecies_;
        scalar grotthussPHReference_;
        scalar grotthussPHWidth_;
        scalar grotthussPHMinFactor_;
        bool grotthussIonicStrengthScaling_;
        wordList grotthussIonicStrengthSpecies_;
        scalarList grotthussIonicStrengthCharges_;
        scalar grotthussIonicStrengthReference_;
        scalar grotthussIonicStrengthExponent_;
        scalar grotthussIonicStrengthMin_;
        scalar grotthussIonicStrengthMinFactor_;

        //Total diffusion flux
        surfaceScalarField phiD_;

        //mass interfacial flux
        volScalarField Mflux_;

public:

    //- Runtime type information
    TypeName("basicTwoPhaseMultiComponentTransportMixture");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicTwoPhaseMultiComponentTransportMixture,
        fvMesh,
        (const fvMesh& mesh,
         const volScalarField& alpha1
        ),
        (mesh, alpha1)
    );


    // Constructors

        //- Construct from mesh and object registry
        basicTwoPhaseMultiComponentTransportMixture
        (
            const fvMesh&,
            const volScalarField& alpha1
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<basicTwoPhaseMultiComponentTransportMixture> New(const fvMesh&, const volScalarField& alpha1);


    //- Destructor
    virtual ~basicTwoPhaseMultiComponentTransportMixture();


    // Member functions

        // Access to phase parameters

        scalar cYi() const
        {
            return cYi_;
        }
        //- Return the species molecular diffusion list for phase 1
        const PtrList<dimensionedScalar>& D1Y() const
        {
            return D1Y_;
        }

        //- Return the species molecular diffusion for phase 1
        dimensionedScalar D1Y(const label i) const
        {
            return D1Y_[i];
        }

        //- Return the species molecular diffusion list for phase 2
        const PtrList<dimensionedScalar>& D2Y() const
        {
            return D2Y_;
        }

        //- Return the species molecular diffusion for phase 2
        dimensionedScalar D2Y(const label i) const
        {
            return D2Y_[i];
        }

        //- Return the species Henry's constant
        dimensionedScalar HY(const label i) const
        {
            return HY_[i];
        }

        //- Return molecular weight
        dimensionedScalar Mw(const label i) const
        {
            return Mw_[i];
        }

        //- Return mass interfacial flux
        const volScalarField Mflux() const
        {
            return Mflux_;
        }

        //- Return mass interfacial flux
        const surfaceScalarField phiD() const
        {
            return phiD_;
        }


        // calculate and return mean diffusion coefficient
        surfaceScalarField DmY(const label i) const
        {
            if
            (
                grotthussActive_[i]
             && (
                    grotthussSaturationSuppression_
                 || grotthussPHScaling_
                 || grotthussIonicStrengthScaling_
                )
            )
            {
                volScalarField structuralFactor
                (
                    IOobject
                    (
                        "grotthussStructuralFactor",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("one", dimless, 1.0)
                );

                if (grotthussSaturationSuppression_)
                {
                    volScalarField alphaEff
                    (
                        min
                        (
                            max
                            (
                                (alpha1_ - dimensionedScalar("minWaterSaturation", dimless, grotthussMinWaterSaturation_))
                               /(dimensionedScalar("oneMinusMinWaterSaturation", dimless, 1.0 - grotthussMinWaterSaturation_)),
                                dimensionedScalar("zero", dimless, 0.0)
                            ),
                            dimensionedScalar("one", dimless, 1.0)
                        )
                    );

                    structuralFactor *= pow(alphaEff, grotthussSaturationExponent_);
                }

                if (grotthussPHScaling_)
                {
                    label hplusIndex(-1);
                    forAll(species_, speciesI)
                    {
                        if (species_[speciesI] == grotthussHplusSpecies_)
                        {
                            hplusIndex = speciesI;
                        }
                    }

                    if (hplusIndex < 0)
                    {
                        FatalErrorInFunction
                            << "grotthussTransport pH scaling requested but species "
                            << grotthussHplusSpecies_ << " was not found."
                            << exit(FatalError);
                    }

                    volScalarField pHDistance
                    (
                        mag
                        (
                           -log10
                            (
                                max
                                (
                                    Y_[hplusIndex],
                                    dimensionedScalar
                                    (
                                        "hplusFloor",
                                        Y_[hplusIndex].dimensions(),
                                        1e-30
                                    )
                                )
                               /dimensionedScalar
                                (
                                    "oneMolar",
                                    Y_[hplusIndex].dimensions(),
                                    1.0
                                )
                            )
                          - dimensionedScalar("pHReference", dimless, grotthussPHReference_)
                        )
                    );

                    structuralFactor *=
                        min
                        (
                            dimensionedScalar("one", dimless, 1.0),
                            max
                            (
                                dimensionedScalar("pHMinFactor", dimless, grotthussPHMinFactor_),
                                exp
                                (
                                   -pHDistance
                                   /dimensionedScalar("pHWidth", dimless, grotthussPHWidth_)
                                )
                            )
                        );
                }

                if (grotthussIonicStrengthScaling_)
                {
                    volScalarField ionicStrength
                    (
                        IOobject
                        (
                            "grotthussIonicStrength",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar("zero", Y_[i].dimensions(), 0.0)
                    );

                    forAll(grotthussIonicStrengthSpecies_, ionicI)
                    {
                        label ionicSpeciesIndex(-1);
                        forAll(species_, speciesI)
                        {
                            if (species_[speciesI] == grotthussIonicStrengthSpecies_[ionicI])
                            {
                                ionicSpeciesIndex = speciesI;
                            }
                        }

                        if (ionicSpeciesIndex < 0)
                        {
                            FatalErrorInFunction
                                << "grotthussTransport ionic strength species "
                                << grotthussIonicStrengthSpecies_[ionicI]
                                << " was not found."
                                << exit(FatalError);
                        }

                        ionicStrength +=
                            0.5
                           *sqr(grotthussIonicStrengthCharges_[ionicI])
                           *max
                            (
                                Y_[ionicSpeciesIndex],
                                dimensionedScalar
                                (
                                    "zero",
                                    Y_[ionicSpeciesIndex].dimensions(),
                                    0.0
                                )
                            );
                    }

                    structuralFactor *=
                        min
                        (
                            dimensionedScalar("one", dimless, 1.0),
                            max
                            (
                                dimensionedScalar("ionicMinFactor", dimless, grotthussIonicStrengthMinFactor_),
                                pow
                                (
                                    dimensionedScalar
                                    (
                                        "ionicReference",
                                        ionicStrength.dimensions(),
                                        grotthussIonicStrengthReference_
                                    )
                                   /max
                                    (
                                        ionicStrength,
                                        dimensionedScalar
                                        (
                                            "ionicMin",
                                            ionicStrength.dimensions(),
                                            grotthussIonicStrengthMin_
                                        )
                                    ),
                                    grotthussIonicStrengthExponent_
                                )
                            )
                        );
                }

                volScalarField D1Eff
                (
                    dimensionedScalar("Dvehicle", D1Y_[i].dimensions(), grotthussDvehicle_[i])
                  + dimensionedScalar("Dstructural", D1Y_[i].dimensions(), grotthussDstructural_[i])*structuralFactor
                );

                return fvc::interpolate(D1Eff*alpha1_+HY_[i]*D2Y_[i]*(1-alpha1_))/fvc::interpolate(alpha1_+HY_[i]*(1-alpha1_));
            }

            return fvc::interpolate(D1Y_[i]*alpha1_+HY_[i]*D2Y_[i]*(1-alpha1_))/fvc::interpolate(alpha1_+HY_[i]*(1-alpha1_));
        } 

        //calculate and return the henry transfer flux
        surfaceScalarField phiH(const label i) const;

        //calculate and return the upwind henry transfer flux
        surfaceScalarField phiHUp(const label i) const;

        //calculate and return the downwind henry transfer flux
        surfaceScalarField phiHDown(const label i) const;

        //calculate compression coefficient
        surfaceScalarField compressionCoeff(const label i);

        //- Update properties
        virtual void correct()=0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
