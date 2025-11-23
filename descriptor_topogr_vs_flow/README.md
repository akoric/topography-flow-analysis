# Topography vs Flow Descriptor Scripts (Paper Companion)

This folder has a small set of Python scripts used to compute and visualize the link between surface/topography descriptors and flow-field descriptors for the paper.

The basic workflow is:

1. Run the **Kendall correlation** analysis (per topography descriptor and flow file).
2. Summarize those results into **heatmaps/rankings** for the paper.
3. Use the **scatter plot** tools for focused visual checks.
4. Optionally generate a **legend** figure used in the paper plots.

Below is a quick summary of each script.

---

## 1. `kendall_correlation_analysis.py`

Single‑descriptor Kendall rank correlation.

- Reads `data/topography_results.csv` and all `data/*.csv` flow files.
- You type a **topography descriptor name** (e.g. `Lx`, `Ly`, `Ra_sw`, etc.).
- For that descriptor it:
  - Orders DEMs by the topography value.
  - For each flow file, grabs level 0 and level 1 values and orders DEMs.
  - Computes Kendall τ between topo ordering and flow ordering at each level.
  - Saves one `results/kendall_report.csv` with τ values and DEM counts.

You normally don’t call this directly; `batch_kendall_analysis.py` drives it.

---

## 2. `batch_kendall_analysis.py`

Batch wrapper around the Kendall analysis.

- Reads `data/topography_results.csv` to get all available **topography descriptors**.
- Lets you choose:
  - `1` → run for **all** descriptors.
  - `2` → pick a subset by number.
  - `3` → just list descriptors and quit.
- For each chosen descriptor it runs `kendall_correlation_analysis.py` once.
- After each run it renames the output to:
  - `results/kendall_report_<descriptor>.csv`.

This is the main entry point to regenerate all Kendall correlation tables used by the paper.

---

## 3. `visualize_kendall_results.py`

Takes the per‑descriptor Kendall reports and builds compact figures and tables.

- Looks for `results/kendall_report_*.csv` produced by the batch script.
- Combines them into one DataFrame.
- Drops fractal‑dimension style descriptors (`FracD_*`) from the main heatmaps.
- Computes:
  - Average |τ| per descriptor (both levels).
  - A simple ranking of descriptors.
  - The best descriptor per flow file and level.
- Outputs to `paper_results/`:
  - `viz_heatmap_level_0.png` (τ, level 0).
  - `viz_heatmap_level_1.png` (τ, level 1).
  - `viz_descriptor_rankings.csv`.
  - `viz_best_descriptors_by_parameter.csv`.

This is what you run when you want the summary plots/tables that appear in the paper.

---

## 4. `descriptor_scatter_plots_interactive.py`

Interactive scatter plots for **one** topo descriptor vs **one** flow descriptor across all directions.

- Uses `data/topography_results.csv` and `data/U_*.csv`, `V_*.csv`, `W_*.csv` flow files.
- You pick:
  - A topography descriptor.
  - A flow parameter type (e.g. `mean`, `kurt`, `skew`, `std`, `var`).
- For all directions U, V, W and for levels 0 and 1, it makes scatter plots of:
  - x: topo descriptor value.
  - y: flow descriptor value.
- Saves into `paper_results/<descriptor>_vs_<flow_param>/` as PNGs.

Use this when you want to visually sanity‑check the Kendall results for a particular descriptor/flow combo.

---

## 5. `descriptor_scatter_plots_sw.py`

Scatter plots with **error bars** and multiple flow parameters.

- Still uses `data/topography_results.csv` plus multiple `data/<flow_param>_results.csv` files.
- Also reads roughness CSVs like `data/roughness_U_mean.csv` to estimate error bars.
- For a chosen descriptor and a list of flow params (e.g. `['U_mean','V_mean','W_mean']`) it:
  - Plots level 0 and level 1 separately.
  - Uses markers for DEMs and vertical error bars from roughness.
- Saves PNGs (level 0 and 1) in the path you give via `save_path`.

This is mainly for the paper’s more detailed scatter figures that show variability.

---

## 6. `generate_legend.py`

Utility to generate a small legend figure that matches marker/color encoding used in the scatter plots.

- Creates a compact figure with:
  - Marker shapes → DEM IDs.
  - Colors → flow directions (U, V, W).
- Saves `paper_results/legend_compact.png`.

You run this once to regenerate the legend art used in the paper figures.

---

## Minimal usage recap

From this folder:

```bash
# 1) Run all Kendall correlations and write per‑descriptor CSVs
python batch_kendall_analysis.py

# 2) Build heatmaps + ranking tables for the paper
python visualize_kendall_results.py

# 3) Optional: interactive scatter plots for one descriptor vs one flow type
python descriptor_scatter_plots_interactive.py

# 4) Optional: scatter plots with error bars for selected flow params
# (typically called from a small helper script or a notebook)

# 5) Optional: regenerate legend figure
python generate_legend.py
```

All scripts assume the `data/` and `results/` folders have already been populated from the upstream simulation/measurement pipeline used in the paper.
