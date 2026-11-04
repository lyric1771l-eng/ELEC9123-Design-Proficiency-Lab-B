# ELEC9123 Design Proficiency - Lab B (2025)

**Wireless Channel Performance Simulation**

MATLAB study of BPSK bit error rate (BER) and outage probability under noise, fading and distance-dependent path loss. Prepared by Yan Liang for ELEC9123 Design Proficiency (Task B), UNSW, 2025.

## Scope

- Monte Carlo BPSK transmission, detection and error counting.
- AWGN and Laplacian noise; Laplacian PDF/CDF checks.
- Rayleigh and Rician fading with theoretical and simulated BER/outage curves.
- Random user placement with path loss and per-user/average performance plots.
- Numerical integration and the Marcum Q-function for Rician performance analysis.

## Entry point

`BL5.m` is the author-confirmed final version. Earlier iterations are intentionally not included. Its three helper functions are defined within the script, so no other coursework `.m` files are needed.

### Run

Open this directory in MATLAB and run:

```matlab
BL5
```

Use MATLAB R2024a or newer for local functions interleaved with script statements. The script uses `marcumq`; ensure the toolbox providing it is installed (`which marcumq`). The MATLAB version used for the original coursework was not recorded here.

The script clears the workspace and closes figures. Initial settings include 1,000,000 bits and an Eb/N0 sweep from 0 to 15 dB. The final deployment section currently sets `U=1`; it should not be described as a large-scale multi-user experiment. Computation and numerical integrations may take time.

## Outputs and interpretation

The script produces BER comparisons, noise-distribution checks, outage curves, and user-deployment plots. These are simulated link-level results, not field measurements or RF hardware tests. Coherent fading detection divides by the channel coefficient, assuming channel knowledge at the receiver despite an inconsistent earlier comment mentioning no CSI.

## Reproducibility and limitations

- Original MATLAB code is preserved byte-for-byte; see `SOURCE_MANIFEST.json`.
- No MATLAB execution or new numerical validation was performed during GitHub packaging. No new performance claims are made.
- The script does not set a random seed. Monte Carlo curves vary between runs; a measured zero BER with finite samples is not proof of zero error probability.
- Review noise normalization across models and numerical integration behavior before using results for engineering decisions.
- Original reports and student-identifying cover pages are not included. This repository presents the final code with usage documentation, rather than every development draft.

## Skills demonstrated

MATLAB · Wireless communications · BPSK · Monte Carlo simulation · Fading channels · BER/outage analysis
