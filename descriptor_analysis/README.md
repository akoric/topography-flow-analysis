# Descriptor Analysis Tools

Scripts for analyzing surface roughness descriptors from DEM velocity field data.

## How It Works

The programs analyze how surface roughness descriptors change over time for different DEM samples. Each descriptor (like Ra, Rq, fractal dimension, etc.) captures different aspects of the velocity field roughness.

**Basic Analysis**: For each descriptor, the code plots values across frames for all DEM samples. This shows temporal evolution and lets you visually compare how different materials behave.

**Slope-Intercept Method**: The advanced analysis fits a linear trend (slope/intercept) for each DEM's descriptor values over time. Then it calculates:
- Separation score: how far apart the DEMs are in slope-intercept space (distinctness)
- Consistency score: how tight the trends are (low variance = high consistency)

Descriptors with high separation and consistency are the best at distinguishing between different materials. The ranking helps identify which descriptors are most useful for characterizing DEM behavior.

## Files

**analyze_descriptors.py** - Interactive tool to analyze a single CSV file. Pick which descriptors to visualize and it generates scatter plots showing how each descriptor varies across DEM samples.

**batch_analyze_descriptors.py** - Automatically processes all CSV files in the data/ directory. Useful for generating plots for multiple velocity components at once.

**consistency_distinctness_analysis.py** - Runs slope-intercept analysis to rank descriptors by how well they separate different DEMs. Creates ranked plots and saves results to CSV.

## Usage

Put your results CSV files in the `data/` directory. Files should have columns for DEM, frame, and descriptor values.

Run any script:
```
python analyze_descriptors.py
python batch_analyze_descriptors.py
python consistency_distinctness_analysis.py
```

Plots get saved to `plots_*` directories.

## Data Format

CSV files need these columns:
- DEM (sample identifier)
- frame (time index)
- descriptor columns (Ra, Rq, Sk, Ku, FracD, Hurst, Lx, Ly, Lk, etc.)

## Dependencies

numpy, pandas, matplotlib, seaborn, scipy, statsmodels
