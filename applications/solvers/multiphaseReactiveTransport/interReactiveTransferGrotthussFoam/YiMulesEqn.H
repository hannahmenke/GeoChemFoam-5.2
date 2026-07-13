{
    word alpharScheme("div(phirb,alpha)");
    word YiScheme("div(phi,Yi)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(speciesMixture.cYi()*mag(phi/mesh.magSf()));

    surfaceScalarField phir(phic*mixture.nHatf());

    basicTwoPhaseMultiComponentTransportMixture& speciesMixture =pSpeciesMixture();
    const speciesTable& solutionSpecies = speciesMixture.species();

    surfaceScalarField phiComp = fvc::flux
    (
        -fvc::flux(-phir, alpha2, alpharScheme),
        alpha1,
        alpharScheme
    );

    autoPtr<volScalarField> h2ActivationH2AqueousSnapshot;
    autoPtr<volScalarField> h2ActivationWeight;
    autoPtr<volScalarField> h2ActivationVolumetricRate;
    autoPtr<volScalarField> h2ActivationExtent;

    if (h2ActivationEnabled)
    {
        h2ActivationH2AqueousSnapshot.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "h2ActivationH2AqueousSnapshot",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                speciesMixture.Y(h2ActivationH2Index)
               /(
                    alpha1
                  + speciesMixture.HY(h2ActivationH2Index)
                   *(dimensionedScalar("one", dimless, 1.0) - alpha1)
                )
            )
        );

        volScalarField h2ActivationSourceWeight
        (
            IOobject
            (
                "h2ActivationSourceWeight",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        volScalarField h2ActivationRateField
        (
            IOobject
            (
                "h2ActivationRateFieldScratch",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        );

        if (h2ActivationMode == "vofInterface")
        {
            volScalarField gradAlphaMag(mag(fvc::grad(alpha1)));
            dimensionedScalar gradAlphaScale
            (
                "h2ActivationGradAlphaScale",
                gradAlphaMag.dimensions(),
                max(gMax(gradAlphaMag), SMALL)
            );

            h2ActivationSourceWeight =
                pos(alpha1 - h2ActivationMinAlpha)
               *pos(h2ActivationMaxAlpha - alpha1)
               *alpha1
               *gradAlphaMag/gradAlphaScale;
            h2ActivationRateField =
                h2ActivationRate*h2ActivationSourceWeight;
        }
        else if (h2ActivationMode == "reactiveWall")
        {
            volScalarField h2ActivationWallAreaPerVolume
            (
                IOobject
                (
                    "h2ActivationWallAreaPerVolume",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimLength, 0.0)
            );

            const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
            forAll(h2ActivationReactivePatches, patchI)
            {
                const word& patchName = h2ActivationReactivePatches[patchI];
                const label patchID = boundaryMesh.findPatchID(patchName);

                if (patchID < 0)
                {
                    FatalErrorInFunction
                        << "h2Activation reactive patch " << patchName
                        << " was not found in the mesh boundary."
                        << abort(FatalError);
                }

                const labelUList& faceCells = boundaryMesh[patchID].faceCells();
                const scalarField& faceAreas = mesh.magSf().boundaryField()[patchID];
                forAll(faceCells, faceI)
                {
                    h2ActivationSourceWeight[faceCells[faceI]] = 1.0;
                    h2ActivationWallAreaPerVolume[faceCells[faceI]] +=
                        faceAreas[faceI]/mesh.V()[faceCells[faceI]];
                }
            }

            h2ActivationSourceWeight *=
                pos(alpha1 - h2ActivationMinAlpha)*alpha1;
            h2ActivationSourceWeight.correctBoundaryConditions();

            if (h2ActivationUseReactiveWallAreaScaling)
            {
                h2ActivationRateField =
                    h2ActivationSurfaceRate
                   *h2ActivationWallAreaPerVolume
                   *h2ActivationSourceWeight;
            }
            else
            {
                h2ActivationRateField =
                    h2ActivationRate*h2ActivationSourceWeight;
            }
        }

        h2ActivationWeight.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "h2ActivationWeight",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                h2ActivationSourceWeight
            )
        );
        h2ActivationVolumetricRate.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "h2ActivationVolumetricRate",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                h2ActivationRateField
            )
        );

        h2ActivationExtent.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "h2ActivationExtent",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                h2ActivationVolumetricRate()
               *h2ActivationH2AqueousSnapshot()
            )
        );

        volScalarField h2ActivationExtentLimit
        (
            IOobject
            (
                "h2ActivationExtentLimit",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            max
            (
                h2ActivationH2AqueousSnapshot(),
                dimensionedScalar("zero", h2ActivationH2AqueousSnapshot().dimensions(), 0.0)
            )/(max(-h2ActivationH2Coefficient, SMALL)*runTime.deltaT())
        );

        if (h2ActivationProtonYield < 0.0)
        {
            const volScalarField& limitingYi =
                speciesMixture.Y(h2ActivationProtonIndex);
            h2ActivationExtentLimit = min
            (
                h2ActivationExtentLimit,
                max
                (
                    limitingYi,
                    dimensionedScalar("zero", limitingYi.dimensions(), 0.0)
                )/(-h2ActivationProtonYield*runTime.deltaT())
            );
        }

        if (h2ActivationHydroxideYield < 0.0)
        {
            const volScalarField& limitingYi =
                speciesMixture.Y(h2ActivationHydroxideIndex);
            h2ActivationExtentLimit = min
            (
                h2ActivationExtentLimit,
                max
                (
                    limitingYi,
                    dimensionedScalar("zero", limitingYi.dimensions(), 0.0)
                )/(-h2ActivationHydroxideYield*runTime.deltaT())
            );
        }

        forAll(h2ActivationYieldSpecies, yieldI)
        {
            const scalar yieldCoeff = h2ActivationYieldCoefficients[yieldI];
            if (yieldCoeff < 0.0)
            {
                const volScalarField& limitingYi =
                    speciesMixture.Y(h2ActivationYieldSpeciesIndices[yieldI]);
                h2ActivationExtentLimit = min
                (
                    h2ActivationExtentLimit,
                    max
                    (
                        limitingYi,
                        dimensionedScalar("zero", limitingYi.dimensions(), 0.0)
                    )/(-yieldCoeff*runTime.deltaT())
                );
            }
        }

        h2ActivationExtent() = min(h2ActivationExtent(), h2ActivationExtentLimit);
        h2ActivationExtent().max
        (
            dimensionedScalar("zero", h2ActivationExtent().dimensions(), 0.0)
        );
        h2ActivationExtent().correctBoundaryConditions();
    }

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);

        //compute diffusion coefficient
        surfaceScalarField D = speciesMixture.DmY(i);

        //calculate Henry's transfer flux
        surfaceScalarField phiH = speciesMixture.phiH(i);
        surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
        surfaceScalarField phiHDown = speciesMixture.phiHDown(i);

        scalar maxYi = max(gMax(Yi.oldTime()),gMax(Yi.oldTime().boundaryField()))+1e-30;

        Yi.oldTime() == Yi.oldTime()/maxYi;
        Yi == Yi/maxYi;

        surfaceScalarField phiDY = -D*fvc::snGrad(Yi)*mesh.magSf();

        if (phiHScheme == "Gauss linear") phiDY += fvc::flux(phiH,Yi,"div(phiH,Yi)");
        else if (phiHScheme == "Gauss upwind")
        {
            phiDY += fvc::flux(phiHUp,Yi,"div(phiH,Yi)")
            + fvc::flux(phiHDown,Yi,"div(phiH,Yi)");
        }
        else
        {
            Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
                            << endl
                            << abort(FatalError);
        }

        volScalarField SpY = -fvc::div(phiDY);

        if (h2ActivationEnabled)
        {
            scalar activationCoeff = 0.0;
            if (solutionSpecies[i] == h2ActivationH2Species)
            {
                activationCoeff = h2ActivationH2Coefficient;
            }
            else if (solutionSpecies[i] == h2ActivationProtonSpecies)
            {
                activationCoeff = h2ActivationProtonYield;
            }
            else if (solutionSpecies[i] == h2ActivationHydroxideSpecies)
            {
                activationCoeff = h2ActivationHydroxideYield;
            }
            forAll(h2ActivationYieldSpecies, yieldI)
            {
                if (solutionSpecies[i] == h2ActivationYieldSpecies[yieldI])
                {
                    activationCoeff += h2ActivationYieldCoefficients[yieldI];
                }
            }

            if (activationCoeff != 0.0)
            {
                SpY +=
                    activationCoeff
                   *h2ActivationExtent()/maxYi;
            }
        }

        tmp<surfaceScalarField> tYiPhi1Un
        (
            fvc::flux
            (
                phi,
                Yi,
                YiScheme
            )
            +   phiComp*speciesMixture.compressionCoeff(i)
        );

        {
            surfaceScalarField YiPhi10 = tYiPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                Yi,
                phi,
                YiPhi10,
                zeroField(),
                SpY,
                oneField(),
                zeroField()
            );
        }

        Yi.oldTime() == Yi.oldTime()*maxYi;
        Yi == Yi*maxYi;

        Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
        << "  Min(Yi) = " << gMin(Yi.internalField())
        << "  Max(Yi) = " << gMax(Yi.internalField())
        << endl;

    }
 }
