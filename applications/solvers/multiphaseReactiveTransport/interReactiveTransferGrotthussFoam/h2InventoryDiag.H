{
    const speciesTable& h2DiagSpecies = speciesMixture.species();
    label h2DiagIndex = -1;

    forAll(h2DiagSpecies, h2DiagI)
    {
        if (h2DiagSpecies[h2DiagI] == "H2")
        {
            h2DiagIndex = h2DiagI;
            break;
        }
    }

    if (h2DiagIndex >= 0)
    {
        const volScalarField& h2DiagField = speciesMixture.Y(h2DiagIndex);
        const scalar h2DiagVolume = gSum(mesh.V());
        const scalar h2DiagInventory =
            h2DiagField.weightedAverage(mesh.V()).value()*h2DiagVolume;

        Info<< "H2InventoryDiag " << h2InventoryDiagStage
            << " time " << runTime.timeName()
            << " inventory_molL_m3 " << h2DiagInventory
            << " avg_molL " << h2DiagField.weightedAverage(mesh.V()).value()
            << " min_molL " << gMin(h2DiagField.internalField())
            << " max_molL " << gMax(h2DiagField.internalField())
            << nl;
    }
}
