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

Application
    interReactiveTransferGrotthussFoam

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "phreeqcMixture.H"
#include "inertMultiComponentMixture.H"
#include "basicTwoPhaseMultiComponentTransportMixture.H"
#include "twoPhaseMultiComponentTransferMixture.H"
#include "fvOptions.H"
#include <CorrectPhi.H>
#include "fvcSmooth.H"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//define monitor time indexes
#define ifMonitor  if( runTime.timeIndex()%10== 0 )

static Foam::scalar parsedSpeciesCharge
(
    const Foam::word& species,
    const bool fatalOnUnparseable = true
)
{
    const std::string speciesName(species.c_str());
    const std::string::size_type plusPos = speciesName.rfind('+');
    const std::string::size_type minusPos = speciesName.rfind('-');

    if (plusPos == std::string::npos && minusPos == std::string::npos)
    {
        return 0.0;
    }

    const bool hasPlus = plusPos != std::string::npos;
    const bool hasMinus = minusPos != std::string::npos;
    const std::string::size_type signPos =
        hasPlus && (!hasMinus || plusPos > minusPos) ? plusPos : minusPos;
    const Foam::scalar sign = speciesName[signPos] == '+' ? 1.0 : -1.0;
    const std::string suffix = speciesName.substr(signPos + 1);
    Foam::label magnitude = 1;

    if (!suffix.empty())
    {
        magnitude = 0;
        for (std::string::size_type charI = 0; charI < suffix.size(); ++charI)
        {
            if (!std::isdigit(static_cast<unsigned char>(suffix[charI])))
            {
                if (fatalOnUnparseable)
                {
                    FatalErrorInFunction
                        << "Cannot parse charge suffix for species " << species
                        << ". Use names like H+, Red-, Fe+3, or CO3-2, "
                        << "or define a neutral species without plus or "
                        << "minus characters."
                        << abort(FatalError);
                }
                return 0.0;
            }
            magnitude = 10*magnitude + suffix[charI] - '0';
        }

        if (magnitude == 0)
        {
            if (fatalOnUnparseable)
            {
                FatalErrorInFunction
                    << "Cannot parse zero charge magnitude for species "
                    << species
                    << abort(FatalError);
            }
            return 0.0;
        }
    }

    return sign*magnitude;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing. \n"
        "with species transport using CST method, phreeqc reactive mixture,"
        " and pH output from H+"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"

            #include "YiMulesEqn.H"
            #include "alphaEqnSubCycle.H"

            gradalpha1 = mag(fvc::grad(alpha1));

            mixture.correct();
            if (h2InventoryDiagnostics)
            {
                const char* h2InventoryDiagStage = "preChem";
                #include "h2InventoryDiag.H"
            }
            #define chargeField chargePreChem
            #include "chargeEqn.H"
            #undef chargeField
            speciesMixture.correct();
            if (h2InventoryDiagnostics)
            {
                const char* h2InventoryDiagStage = "postChem";
                #include "h2InventoryDiag.H"
            }
            #define chargeField chargePostChem
            #include "chargeEqn.H"
            #undef chargeField
            #include "pHEqn.H"

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        //monitor average and max velocity
        ifMonitor
        {
                        Info << "\n         Umax = " << max(mag(U)).value() << " m/s  "
                        << "Uavg = " << mag(average(U)).value() << " m/s";
        }


        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
