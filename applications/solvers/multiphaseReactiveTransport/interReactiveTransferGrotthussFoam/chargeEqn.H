{
    chargeField = dimensionedScalar
    (
        "zeroCharge",
        chargeField.dimensions(),
        0.0
    );

    const speciesTable& chargeSpecies = speciesMixture.species();

    forAll(chargeSpecies, i)
    {
        scalar z = parsedSpeciesCharge(chargeSpecies[i]);

        if (z != 0.0)
        {
            chargeField += z*speciesMixture.Y(i);
        }
    }
}
