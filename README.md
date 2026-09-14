# GARCH_kriging

Reproducible R and Quarto code for heterogeneous spatiotemporal GARCH models, with a companion tutorial on spatial volatility prediction.

The repository accompanies:

**Atika Aouri and Philipp Otto (2025).** *A Heterogeneous Spatiotemporal GARCH Model: A Predictive Framework for Volatility in Financial Networks.* [arXiv:2508.20101](https://arxiv.org/abs/2508.20101).

## Contents

- Local estimation of spatially varying GARCH parameters.
- Spatial prediction using squared-innovation best linear prediction and Gaussian conditional moments.
- Simulation examples and applications to air pollution and financial data.
- Quarto reproducibility materials and a worked tutorial.

## Getting started

1. Clone or download the repository, preserving its directory structure.
2. Install R, Quarto, and the R packages required by the scripts.
3. Follow the data-loading instructions in the relevant Quarto document and ensure that its input files are available.
4. Run the required computation blocks, then render the document to produce the figures and tables.

**Rendering is not the same as rerunning every analysis.** Computationally intensive blocks marked `eval: false` must be run manually to regenerate their saved results. Subsequent blocks load those results when the document is rendered. Fixed random seeds are specified in the code; some estimation routines support parallel execution.

The reduced simulation is an illustrative example, not a substitute for the full Monte Carlo study. The tutorial explains the workflow, while the manuscript reproducibility materials document the corresponding analyses.

## Citation

If you use the software or methodology in your research, please cite:

> Aouri, A., & Otto, P. (2025). A Heterogeneous Spatiotemporal GARCH Model: A Predictive Framework for Volatility in Financial Networks. arXiv. https://doi.org/10.48550/arXiv.2508.20101

Citation metadata are available in [CITATION.cff](CITATION.cff). Please also record the repository version or commit used for your analysis.

## Questions and contributions

Questions, suggestions, and reproducibility reports are welcome through the issue tracker. Please include the relevant script, software versions, and any error message.

Third-party data and materials remain subject to their original licenses and access conditions.
