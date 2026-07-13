# Grotthuss GeoChemFoam Handover

## Purpose

This branch adds a Grotthuss-aware reactive-transfer solver to GeoChemFoam 5.2
for pore-scale hydrogen dissolution and reactive-transport studies. It is the
solver dependency for the companion `GrotthussMechanism` project, which holds
the cases, validation tools, execution scripts, reports, and manuscript.

This document records the state on 2026-07-13 and the work still required. It
does not claim that the remaining production cases have been run.

## Repositories

- Solver fork: `https://github.com/hannahmenke/GeoChemFoam-5.2`
- Solver branch: `feature/grotthuss-reactive-transport`
- Upstream base: `GeoChemFoam/GeoChemFoam-5.2`, commit `33c0fb3`
- Companion project: `GrotthussMechanism`, currently a separate local Git
  repository

The solver branch must be pushed only to the personal fork unless an upstream
contribution is explicitly authorized later.

## Implemented Solver Work

The branch contains:

- `interReactiveTransferGrotthussFoam`, forked from the multiphase reactive
  transfer solver;
- separate vehicle and structural contributions to proton and hydroxide
  mobility;
- pH and activity-based `pHActivity` output;
- charge diagnostics before and after PHREEQC coupling;
- H2 inventory diagnostics;
- bounded H2 activation source modes, including VOF-interface and reactive-wall
  activation;
- acceptor and acid/base availability limiting;
- sulfate-analogue and PHREEQC-native sulfate-to-sulfide field support in the
  companion cases;
- build integration through
  `applications/solvers/multiphaseReactiveTransport/Allwmake`.

The implementation touches the new solver plus:

- `basicTwoPhaseMultiComponentTransportMixture.C/.H`;
- `phreeqcMixture.C/.H`.

## Build And Smoke Check

The companion project provides `build_grotthuss_in_docker.sh`. From its root:

```bash
./build_grotthuss_in_docker.sh geochemfoam/geochemfoam-5.2:latest
```

The script mounts this source checkout read-only, copies the modified sources
into the image workspace, rebuilds `reactionThermo`, and builds
`interReactiveTransferGrotthussFoam`.

On a new execution host, first confirm:

1. the Docker image is available;
2. `wmake libso src/thermophysicalModels/reactionThermo` succeeds;
3. `wmake applications/solvers/multiphaseReactiveTransport/interReactiveTransferGrotthussFoam`
   succeeds;
4. the executable resolves in the configured OpenFOAM environment;
5. one short null case reproduces stock H2 dissolution when structural
   transport and activation are disabled.

Do not begin the production matrix until these checks pass. Record the image
digest, source commit, compiler output, and executable checksum for provenance.

## Companion Project Transfer

The companion project Git repository tracks source code, case templates,
scripts, manuscript source, and curated PNG figures. Its `.gitignore`
intentionally excludes solver outputs, PDFs, numerical datasets, raw images,
meshes, and `validation_runs/`.

Consequently, cloning Git repositories alone will not transfer the prepared
production cases or source geometry datasets. Before using the larger machine,
do one of the following:

1. transfer the ignored `validation_runs/maartje_*` directories and required
   raw geometry assets separately, preserving checksums; or
2. transfer the original geometry assets and regenerate all prepared cases with
   the companion staging tools.

After transfer, rerun the preflight tools and compare checksums. Never treat an
incomplete copied directory as a runnable case.

## Current Accepted Evidence

The following work is complete and should not be rerun merely to recreate a
plot:

- focused closed-H2 Henry validation, including phase partition and inventory
  closure;
- stock-solver versus Grotthuss-solver null equivalence for neutral H2
  dissolution with activation disabled;
- 1D analytical diffusion checks and mesh/timestep sensitivity;
- direct full-mobility versus vehicle-only moment and arrival-time comparison;
- Case A electrolyte-support sweep;
- Case B carbonate-buffer sweep;
- Case C Peclet sweep and refinement diagnosis;
- strict body-fitted mesh acceptance for both exact Maartje chips;
- PHREEQC batch and runtime smoke support for sulfate-to-sulfide chemistry.

The homogeneous exact mesh contains 960,508 cells. The cleaned heterogeneous
mesh contains 667,535 cells. Both passed strict
`checkMesh -allTopology -allGeometry` with `Mesh OK.` and have exact cached cell
centres and volumes.

## Prepared But Not Executed

All active simulations were stopped on 2026-07-13. No production run should be
assumed active or complete.

Six passive exact-chip pairs are prepared for later execution:

- homogeneous Pe = 0.3, 1, and 3;
- heterogeneous Pe = 0.3, 1, and 3.

Each pair contains:

- an exact uniform-concentration preservation control;
- a standardized two-inlet mixing case;
- ten checkpoints;
- a 17-digit, MPI-safe finite-volume flux ledger;
- reporting and restart scripts.

The archived `*_ledger6digit_rejected` run failed the strict budget threshold
because output precision accumulated a 1.2161e-3 error. The archived
`*_precision17_partial_stopped` run was intentionally stopped. Neither is
acceptance evidence.

## Required Execution Sequence

Follow `maartje_exact_chip_execution_runbook.md` in the companion project. The
gate order is mandatory.

### Gate 1: Reconfirm Mesh Provenance

- Verify mesh checksums after transfer.
- Re-run strict `checkMesh` only if the mesh or platform transfer is uncertain.
- Confirm exact `cellCenters` and `cellVolumes` counts and summed volumes.

### Gate 2: Passive Pilot

- Run the homogeneous Pe = 1 uniform control first.
- Then run its matched mixing case.
- Require exact constant-field preservation.
- Require per-step finite-volume budget closure within 1e-3.
- Generate the passive report and spatial field figure.

### Gate 3: Passive Matrix

- Run the remaining homogeneous Pe cases serially.
- Run the heterogeneous Pe cases serially only after the pilot gate passes.
- Run one heavy Docker job at a time.
- Set `MpiRanks` for the larger host before decomposition and do not change the
  decomposition width when resuming a case.

### Gate 4: Tier 1 Mobility Comparison

- Generate matched full-mobility and vehicle-only cases only after the
  corresponding passive row passes.
- Keep every setting identical except H+ and OH- mobility.
- Require at least five positive common write times.
- Report activity-based pH when available, pre/post-PHREEQC charge diagnostics,
  second moments, arrival times, and paired spatial fields.

### Gate 5: Tier 2 Buffering

- Run homogeneous unbuffered and weak-buffer design controls first.
- Do not create an experiment-labelled branch until measured alkalinity, total
  inorganic carbon, salinity, and initial pH are available.
- Compare chemical observability, not only imposed diffusivity.

### Gate 6: Tier 3 H2 And Sulfate

- Run matched full/vehicle H2-sulfate co-injection in both exact chips.
- Use real `SO4-2` and `HS-` fields.
- Save OH- or alkalinity, sulfide, sulfate, H2, activity-based pH, cumulative
  reaction extent, and charge diagnostics.
- Verify the 4:1:1:1 stoichiometric ledger.
- Treat the source rate as a bounded design sensitivity until site-specific
  microbial kinetics are calibrated.
- Keep `diffusion_acceleration_factor = 1` in flow-defined Peclet cases.

## Runtime Discipline

Before every long job:

1. estimate wall time from the closest measured sibling case;
2. record the estimate and selected MPI width;
3. set a coarse check interval, normally hourly;
4. work on independent analysis or manuscript tasks while it runs;
5. never repeatedly poll the solver.

Only one heavy Docker job should run at a time. Preserve checkpoints and use the
provided resume scripts after interruption.

## Interpretation Limits

The implemented transport is enhanced independent-Fickian H+ and OH- mobility.
It does not solve the Nernst-Planck electric-potential problem and does not
enforce ambipolar transport.

Therefore:

- quantitative Grotthuss claims should use the swamping-electrolyte trace-ion
  regime;
- low-electrolyte results are upper bounds or limitation demonstrations;
- unsupported HCl-like diffusion should approach the ambipolar HCl coefficient,
  which this implementation does not reproduce;
- the candidate OH- vehicle/structural split is a sensitivity partition, not a
  uniquely measured decomposition;
- bubble shrinkage is a neutral-H2 null observable unless activation couples H2
  to acid/base chemistry;
- the mechanism-facing bubble result is the chemical halo or footprint, not a
  claim that Grotthuss transport directly accelerates neutral H2 diffusion.

## Remaining Model Development

The following improvements remain outside the current validated claim:

- Nernst-Planck or another charge-consistent ambipolar transport formulation for
  unsupported electrolytes;
- calibrated sulfate-reduction kinetics from measured experimental conditions;
- sensitivity bounds for activation kinetics and acceptor availability;
- a production 2D pulsed-source moment fit that avoids continuous-source
  variance ambiguity;
- final matched real-field sulfate full/vehicle visual comparisons on both exact
  chips;
- broader 2D mesh and timestep convergence for the final application metrics.

## Reporting And Publication Tasks

After each accepted gate:

- regenerate machine-readable reports and publication figures from the accepted
  case only;
- update `validation_status.md`, `validation_requirements_audit.md`, and the
  validation document;
- rebuild the manuscript PDF on the execution or writing host, but do not commit
  compiled PDFs or raw solver datasets;
- label verification, validation, and application demonstrations separately;
- retain rejected and partial runs as diagnostics, never as headline evidence;
- archive code, manifests, input data, accepted outputs, and checksums in a
  citable repository before submission.

For every bibliography entry, verify authors, title, journal, year, volume,
issue, pages or article number, and DOI against the publisher, Crossref record,
or paper PDF. DOI resolution alone is insufficient. Any numerical calibration
value must be checked directly in the paper text, table, or figure before use.

## Completion Definition

The study is ready for its final application claim only when:

1. all passive conservation and uniform-field gates pass for both chips;
2. matched Tier 1 fields quantify the mobility effect in both geometries;
3. Tier 2 demonstrates when buffering hides that effect;
4. Tier 3 demonstrates the H2-sulfate chemical footprint with bounded and fully
   documented kinetics;
5. charge, species, and reaction ledgers close within declared tolerances;
6. figures are regenerated from accepted cases and the manuscript claims match
   those results;
7. solver source, case inputs, geometry provenance, and accepted result datasets
   are archived with checksums and a persistent identifier.
