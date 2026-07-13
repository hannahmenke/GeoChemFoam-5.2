{
    const speciesTable& solutionSpecies = speciesMixture.species();
    label hPlusIndex = -1;

    forAll(solutionSpecies, i)
    {
        if (solutionSpecies[i] == "H+")
        {
            hPlusIndex = i;
            break;
        }
    }

    if (hPlusIndex >= 0)
    {
        const volScalarField& Hplus = speciesMixture.Y(hPlusIndex);
        const scalar hFloor = 1e-30;

        forAll(pH, celli)
        {
            if (alpha1[celli] > 1e-3)
            {
                pH[celli] = -std::log10(max(Hplus[celli], hFloor));
            }
            else
            {
                pH[celli] = 30.0;
            }
        }

        forAll(pH.boundaryField(), patchi)
        {
            scalarField& pHpatch = pH.boundaryFieldRef()[patchi];
            const scalarField& Hpatch = Hplus.boundaryField()[patchi];
            const scalarField& alphaPatch = alpha1.boundaryField()[patchi];

            forAll(pHpatch, facei)
            {
                if (alphaPatch[facei] > 1e-3)
                {
                    pHpatch[facei] = -std::log10(max(Hpatch[facei], hFloor));
                }
                else
                {
                    pHpatch[facei] = 30.0;
                }
            }
        }
    }
}
