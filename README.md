# Topography–Flow Analysis Scripts (Amar Koric Contributions)

This repository contains analysis scripts developed by **Amar Koric** as part of the paper:

> **"Quantifying Topography–Flow Relationships: Implications for Turbulent Flow and Roughness Characterization"**

This work was carried out under the supervision of **Dr. Leonardo P. Chamorro** in the **Renewable Energy & Turbulent Environment Group** at the **University of Illinois Urbana-Champaign**.

The code here was used to explore and quantify relationships between digital elevation model (DEM) topography and turbulent flow statistics. Not all scripts and not all results produced from these programs will appear in the final version of the paper, but they reflect the broader set of tests and diagnostics that informed the work.

## Data availability

The MATLAB and Python scripts expect DEM fields and processed descriptor tables that are not included in this repository. The underlying data are not yet public because the paper is still in preparation.

Once the paper has been submitted, the datasets and any additional processing details will be made available on request. Until then, paths and filenames in the code should be treated as examples from the working environment rather than a complete, ready-to-run package.

## MATLAB scripts (root folder)

### `parameter_calculation.m`

Computes a suite of spectral and roughness descriptors from DEM-based velocity fields over a specified range of frames for several DEM cases. For each variable (e.g. `U_mean`, `U_std`, `U_var`, and analogous quantities for V and W) it:

- Loads per-DEM moment files.
- Extracts 2D frames and removes all-NaN rows/columns.
- Computes 2D and 1D power spectral densities (PSDs), radial averages, and associated fractal/Hurst exponents.
- Computes covariance-based length scales and basic tribological roughness measures (Ra, Rq, skewness, kurtosis, etc.).
- Aggregates results into tables and writes them to `.mat` and `.csv` files for later analysis.

These outputs feed into the Python descriptor tools in the subdirectories below.

### `comparisons.m`

Performs post-processing and comparison between flow-derived descriptors and DEM-based roughness descriptors. In particular, it:

- Loads selected descriptor tables for a small set of DEMs and identifies representative frames for each DEM.
- Normalizes quantities (e.g. by maximum mean velocity or domain length) to produce non-dimensional combinations.
- Systematically tests combinations of roughness/flow parameters raised to integer powers, and fits simple linear relationships between flow and DEM combinations.
- Records combinations with high coefficient of determination (R²), and produces simple scatter plots for the best-performing relationships.

This script was used as an initial brute-force exploratory study to scan for compact descriptor combinations that could guide and motivate the more focused analyses and methods developed later in the project.

## Python descriptor tools

### `descriptor_analysis/`

Tools for analyzing how various roughness and spectral descriptors behave across DEMs and time:

- `analyze_descriptors.py`: interactive analysis of a single results CSV, with scatter plots of descriptor values across DEMs and frames.
- `batch_analyze_descriptors.py`: batch driver to run descriptor analysis over multiple CSV files in `descriptor_analysis/data/`.
- `consistency_distinctness_analysis.py`: slope–intercept style analysis used to rank descriptors by how consistently and distinctly they separate DEMs.

These scripts were used to identify which descriptors are most reliable and discriminative; only a subset of the generated figures and rankings are expected to appear in the final paper.

### `descriptor_topogr_vs_flow/`

Scripts focused on relating topographic descriptors to flow-field descriptors using rank correlations and scatter plots:

- `batch_kendall_analysis.py` and `kendall_correlation_analysis.py`: compute Kendall rank correlations between topography descriptors and flow descriptors over multiple files.
- `visualize_kendall_results.py`: aggregates per-descriptor correlation results into heatmaps, rankings, and tables suitable for inclusion in the paper.
- `descriptor_scatter_plots_interactive.py` and `descriptor_scatter_plots_sw.py`: generate various scatter plots (with and without error bars) to visualize the relationships between chosen topography and flow descriptors.
- `generate_legend.py`: utility for producing a compact legend figure used alongside the scatter plots.

As with the rest of the repository, not every script, plot, or intermediate result from this directory will be used in the final manuscript, but they document the analysis space that was explored.
