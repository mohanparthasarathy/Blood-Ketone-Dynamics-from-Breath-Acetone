# Multiscale breath-ketone reconstruction

MATLAB implementation of the metabolic, respiratory transport, and optical models used to investigate reconstruction of blood ketone trajectories from repeated simulated breaths. The project combines numerical verification, patient-specific physics-informed neural networks (PINNs), Fisher-information measurement design, and uncertainty analyses. This is research software and does not diagnose disease.

## Requirements

MATLAB R2024b or newer and Deep Learning Toolbox are required for the PINN routines. Full optimization and response-profile runs are computationally intensive. The code uses `ode15s` for the metabolic model and the finite-difference airway system. Random seeds and optimization options are specified in the corresponding functions.

## Data

`synthetic_data/outputs/cohort/synthetic_cohort.mat` contains the provided NIST-based simulated cohort. `synthetic_data/67-64-1-IR.jdx` is the reference acetone spectrum from the NIST Chemistry WebBook, Coblentz No. 08767: https://webbook.nist.gov/cgi/cbook.cgi?ID=C67641&Index=2&Type=IR-SPEC. The reference spectrum supplies a relative shape; absolute scale, measurement noise, and background spectra are simulation assumptions. No human data are included.

## Workflow

Start MATLAB in the repository root. Modules contain routines with the same names, so run stages in fresh MATLAB sessions to avoid path conflicts.

1. Forward verification and the direction-rule analysis: `run_validation`. To regenerate a cohort, add `synthetic_data` to the MATLAB path and run `run_inverse_study`. Generated files are written to `synthetic_data/outputs/cohort`.
2. Patient-specific reconstruction: `run_reconstruction`. The oracle and informed-prior results are written under `patient_pinn/outputs/`. Each condition fits the supplied test patients with temporal holdouts.
3. Measurement selection: add `classical_fisher` to the path and run `run_fisher_full`. The end-to-end measurement-design comparison is run by adding `fisher_pinn_validation` to the path and running `run_fisher_pinn_full`. The latter adds spectral nuisance terms and retrains for each measurement design.
4. Response-profile analysis: after the informed-prior reconstruction, add `pinn_uncertainty` to the path and run `run_uncertainty_full`. The profile envelopes describe sensitivity to the response parameter; they are not calibrated clinical confidence intervals.
5. Robustness tests: add `final_project_pipeline` to the path and run `run_final_full`. Checkpoints are saved inside that module's output directory. Additional Fisher-design profiles are in `fisher_uncertainty`.

The `run_*_smoke_test` functions provide smaller checks before full fitting. `make_history_figure.m`, `make_trajectory_figure.m`, `make_uncertainty_figure.m`, and `synthetic_data/plot_nist_fisher_channels.m` generate figures from saved outputs. The root figure scripts use explicit paths to the cohort, informed-prior fit, and full response-profile results. Publication figures are in `publication/`.

## Organization

- `synthetic_data/`: cohort generation, solvers, input spectrum, and verification.
- `patient_pinn/`: patient-specific full-measurement reconstruction.
- `classical_fisher/`: Fisher sensitivities and measurement selection.
- `fisher_pinn_validation/`: full-grid and reduced-design comparison.
- `pinn_uncertainty/`, `fisher_uncertainty/`: response profiles and sensitivity bands.
- `final_project_pipeline/`: robustness and exploratory alert calculations.

Some modules include their own versions of training functions. This preserves each experiment's implementation. Outputs are generated locally and may vary with numerical environment and optimization. The modeled alert thresholds are exploratory research choices, not clinical guidance.

## Publication release (1.0.0)

The release preserves completed MATLAB experiments and their saved outputs. For the manuscript figures, use `publication/`. To regenerate the corrected uncertainty figure from CSV data, install `publication/requirements.txt` and run `python publication/figure_source/rebuild_uncertainty.py`. Run `python publication/check_direction.py` to reproduce the matched-observation direction scores (139/240 and 211/240). The original recoverability summary evaluates the plateau rule on a different, test-only denominator; it is not the source of the matched-observation comparison.

### Experimental definitions

- The evaluated sparse design uses eight Fisher-selected channels: 1210, 1226, 1450, 1630, 1720, 1736, 1750, and 1900 cm⁻¹.
- The reconstruction comparison uses a fixed within-breath schedule. Its first two samples are 1.35 and 1.85 seconds; the separate new phase search returns 0.85 and 1.65 seconds. Do not substitute the phase-search output when describing the saved comparison.
- `peak_8` denotes a fixed peak-region comparator (1196, 1210, 1226, 1240, 1706, 1720, 1736, 1750 cm⁻¹), not the strongest eight NIST channels.
- Main manuscript acetone errors are means of patient-specific RMSEs. The bootstrap statistics script pools acetone residuals; these are distinct summaries. BHB and slope RMSE pool held-out observations.
- Robustness conditions change both training seeds and response-prior perturbations. Their results describe combined fitting and measurement variation, not isolated causal effects of noise or sampling.
- Profile bands are sensitivity envelopes. The plotted point curve is the minimum-misfit profile refit, not the main informed-prior fit. Coverage is 14/18 for BHB and 17/18 for slope.

### Release provenance

Numerical solver and training routines and saved numerical outputs are retained from the completed run archive. Release edits concern figure routing, figure labels, publication exports, documentation, and the matched-observation score calculation. `MANIFEST.sha256` records release file checksums. The software is distributed under the MIT license; reference data are credited separately. MATLAB modules retain experiment-specific functions with duplicate names, so use a fresh MATLAB session for each stage as directed above.
