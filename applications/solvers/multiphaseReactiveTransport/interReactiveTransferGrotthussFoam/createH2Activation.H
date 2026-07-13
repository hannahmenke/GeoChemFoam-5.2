bool h2ActivationEnabled(false);
word h2ActivationMode("vofInterface");
word h2ActivationH2Species("H2");
word h2ActivationProtonSpecies("H+");
word h2ActivationHydroxideSpecies("OH-");
dimensionedScalar h2ActivationRate
(
    "rateConstant",
    dimless/dimTime,
    0.0
);
dimensionedScalar h2ActivationSurfaceRate
(
    "surfaceRateConstant",
    dimLength/dimTime,
    0.0
);
scalar h2ActivationProtonYield(0.0);
scalar h2ActivationHydroxideYield(0.0);
scalar h2ActivationH2Coefficient(-1.0);
bool h2ActivationAllowUnbalanced(false);
bool h2ActivationUseReactiveWallAreaScaling(false);
scalar h2ActivationMinAlpha(1e-6);
scalar h2ActivationMaxAlpha(1.0 - 1e-6);
wordList h2ActivationYieldSpecies;
scalarList h2ActivationYieldCoefficients;
labelList h2ActivationYieldSpeciesIndices;
wordList h2ActivationReactivePatches;

if (thermoPhysicalProperties.found("h2Activation"))
{
    const dictionary& h2ActivationDict =
        thermoPhysicalProperties.subDict("h2Activation");

    h2ActivationEnabled =
        h2ActivationDict.lookupOrDefault<bool>("enabled", false);
    h2ActivationMode =
        h2ActivationDict.lookupOrDefault<word>("mode", "vofInterface");
    h2ActivationH2Species =
        h2ActivationDict.lookupOrDefault<word>("h2Species", "H2");
    h2ActivationProtonSpecies =
        h2ActivationDict.lookupOrDefault<word>("protonSpecies", "H+");
    h2ActivationHydroxideSpecies =
        h2ActivationDict.lookupOrDefault<word>("hydroxideSpecies", "OH-");
    h2ActivationProtonYield =
        h2ActivationDict.lookupOrDefault<scalar>("protonYield", 0.0);
    h2ActivationHydroxideYield =
        h2ActivationDict.lookupOrDefault<scalar>("hydroxideYield", 0.0);
    h2ActivationH2Coefficient =
        h2ActivationDict.lookupOrDefault<scalar>("h2Coefficient", -1.0);
    if (h2ActivationDict.found("yieldSpecies"))
    {
        h2ActivationDict.lookup("yieldSpecies") >> h2ActivationYieldSpecies;
    }
    if (h2ActivationDict.found("yieldCoefficients"))
    {
        h2ActivationDict.lookup("yieldCoefficients") >> h2ActivationYieldCoefficients;
    }
    h2ActivationAllowUnbalanced =
        h2ActivationDict.lookupOrDefault<bool>("allowUnbalanced", false);
    h2ActivationUseReactiveWallAreaScaling =
        h2ActivationDict.lookupOrDefault<bool>("useReactiveWallAreaScaling", false);
    h2ActivationMinAlpha =
        h2ActivationDict.lookupOrDefault<scalar>("minWaterSaturation", 1e-6);
    h2ActivationMaxAlpha =
        h2ActivationDict.lookupOrDefault<scalar>("maxWaterSaturation", 1.0 - 1e-6);
    if (h2ActivationDict.found("reactivePatches"))
    {
        h2ActivationDict.lookup("reactivePatches") >> h2ActivationReactivePatches;
    }

    if (h2ActivationDict.found("rateConstant"))
    {
        h2ActivationRate = dimensionedScalar("rateConstant", h2ActivationDict);
    }
    if (h2ActivationDict.found("surfaceRateConstant"))
    {
        h2ActivationSurfaceRate =
            dimensionedScalar("surfaceRateConstant", h2ActivationDict);
    }
}

const speciesTable& h2ActivationSpecies = speciesMixture.species();
label h2ActivationH2Index(-1);
label h2ActivationProtonIndex(-1);
label h2ActivationHydroxideIndex(-1);

forAll(h2ActivationSpecies, h2ActivationI)
{
    if (h2ActivationSpecies[h2ActivationI] == h2ActivationH2Species)
    {
        h2ActivationH2Index = h2ActivationI;
    }
    if (h2ActivationSpecies[h2ActivationI] == h2ActivationProtonSpecies)
    {
        h2ActivationProtonIndex = h2ActivationI;
    }
    if (h2ActivationSpecies[h2ActivationI] == h2ActivationHydroxideSpecies)
    {
        h2ActivationHydroxideIndex = h2ActivationI;
    }
}

if (h2ActivationEnabled)
{
    if (h2ActivationMode != "vofInterface" && h2ActivationMode != "reactiveWall")
    {
        FatalErrorInFunction
            << "Unsupported h2Activation mode " << h2ActivationMode
            << ". Supported modes are vofInterface and reactiveWall."
            << abort(FatalError);
    }

    if (h2ActivationMode == "reactiveWall" && h2ActivationReactivePatches.empty())
    {
        FatalErrorInFunction
            << "h2Activation mode reactiveWall requires reactivePatches."
            << abort(FatalError);
    }

    if (h2ActivationH2Index < 0)
    {
        FatalErrorInFunction
            << "h2Activation h2Species " << h2ActivationH2Species
            << " was not found in solutionSpecies."
            << abort(FatalError);
    }

    if (h2ActivationProtonYield < 0.0 && h2ActivationProtonIndex < 0)
    {
        FatalErrorInFunction
            << "h2Activation protonSpecies " << h2ActivationProtonSpecies
            << " was not found in solutionSpecies, but protonYield is negative."
            << abort(FatalError);
    }

    if (h2ActivationHydroxideYield < 0.0 && h2ActivationHydroxideIndex < 0)
    {
        FatalErrorInFunction
            << "h2Activation hydroxideSpecies " << h2ActivationHydroxideSpecies
            << " was not found in solutionSpecies, but hydroxideYield is negative."
            << abort(FatalError);
    }

    const scalar h2ActivationChargeYield =
        h2ActivationProtonYield - h2ActivationHydroxideYield;

    scalar h2ActivationExtraChargeYield(0.0);

    if
    (
        h2ActivationYieldSpecies.size()
     != h2ActivationYieldCoefficients.size()
    )
    {
        FatalErrorInFunction
            << "h2Activation yieldSpecies and yieldCoefficients must have "
            << "the same length. Found "
            << h2ActivationYieldSpecies.size() << " species and "
            << h2ActivationYieldCoefficients.size() << " coefficients."
            << abort(FatalError);
    }

    h2ActivationYieldSpeciesIndices.setSize(h2ActivationYieldSpecies.size(), -1);

    forAll(h2ActivationYieldSpecies, yieldI)
    {
        forAll(h2ActivationSpecies, speciesI)
        {
            if (h2ActivationSpecies[speciesI] == h2ActivationYieldSpecies[yieldI])
            {
                h2ActivationYieldSpeciesIndices[yieldI] = speciesI;
            }
        }

        if (h2ActivationYieldSpeciesIndices[yieldI] < 0)
        {
            FatalErrorInFunction
                << "h2Activation yield species "
                << h2ActivationYieldSpecies[yieldI]
                << " was not found in solutionSpecies."
                << abort(FatalError);
        }

        h2ActivationExtraChargeYield +=
            h2ActivationYieldCoefficients[yieldI]
           *parsedSpeciesCharge(h2ActivationYieldSpecies[yieldI]);
    }

    if
    (
        mag(h2ActivationChargeYield + h2ActivationExtraChargeYield) > SMALL
     && !h2ActivationAllowUnbalanced
    )
    {
        FatalErrorInFunction
            << "h2Activation has nonzero net charge yield "
            << h2ActivationChargeYield + h2ActivationExtraChargeYield
            << ". Set allowUnbalanced true only for documented "
            << "phenomenological sensitivity runs, or define a "
            << "charge-balanced stoichiometric yield list."
            << abort(FatalError);
    }

    Info<< "H2 activation enabled: mode=" << h2ActivationMode
        << ", rateConstant=" << h2ActivationRate
        << ", surfaceRateConstant=" << h2ActivationSurfaceRate
        << ", h2Coefficient=" << h2ActivationH2Coefficient
        << ", protonYield=" << h2ActivationProtonYield
        << ", hydroxideYield=" << h2ActivationHydroxideYield
        << ", extraYieldSpecies=" << h2ActivationYieldSpecies
        << ", extraYieldCoefficients=" << h2ActivationYieldCoefficients
        << ", extraYieldSpeciesIndices=" << h2ActivationYieldSpeciesIndices
        << ", protonSpeciesIndex=" << h2ActivationProtonIndex
        << ", hydroxideSpeciesIndex=" << h2ActivationHydroxideIndex
        << ", reactivePatches=" << h2ActivationReactivePatches
        << ", useReactiveWallAreaScaling="
        << h2ActivationUseReactiveWallAreaScaling
        << ", allowUnbalanced=" << h2ActivationAllowUnbalanced
        << endl;
}
