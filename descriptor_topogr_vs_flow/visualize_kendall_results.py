import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    import seaborn as sns
    HAS_SEABORN = True
except Exception:
    HAS_SEABORN = False

RESULTS_DIR = "results"
PAPER_RESULTS_DIR = "paper_results"
FIGSIZE_WIDE = (14, 8)
FIGSIZE_TALL = (10, 12)
CMAP = "RdBu_r"
 


def load_all_reports(results_dir=RESULTS_DIR):
    report_files = glob.glob(os.path.join(results_dir, "kendall_report_*.csv"))
    if not report_files:
        raise FileNotFoundError(f"No report files found in {results_dir}")

    all_rows = []
    for path in sorted(report_files):
        desc = os.path.basename(path)[len("kendall_report_"):-4]
        df = pd.read_csv(path)
        df = df.copy()
        df["descriptor"] = desc
        all_rows.append(df)

    combined = pd.concat(all_rows, ignore_index=True)
    return combined


def compute_descriptor_rankings(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for desc, group in df.groupby("descriptor"):
        avg_abs_tau_0 = group["tau_level_0"].abs().mean()
        avg_abs_tau_1 = group["tau_level_1"].abs().mean()
        overall = (avg_abs_tau_0 + avg_abs_tau_1) / 2
        taus_all = pd.concat([group["tau_level_0"], group["tau_level_1"]])
        std_all = float(taus_all.std()) if len(taus_all) > 1 else 0.0
        consistency = 1.0 / (1.0 + std_all)
        strong = int((group["tau_level_0"].abs() >= 0.6).sum() + (group["tau_level_1"].abs() >= 0.6).sum())
        perfect = int((group["tau_level_0"].abs() == 1.0).sum() + (group["tau_level_1"].abs() == 1.0).sum())
        rows.append({
            "descriptor": desc,
            "overall_avg_abs_tau": overall,
            "avg_abs_tau_0": avg_abs_tau_0,
            "avg_abs_tau_1": avg_abs_tau_1,
            "consistency": consistency,
            "strong_correlations": strong,
            "perfect_correlations": perfect,
        })
    rank_df = pd.DataFrame(rows).sort_values("overall_avg_abs_tau", ascending=False)
    return rank_df


def compute_best_descriptor_per_param(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for param, sub in df.groupby("file"):
        # Level 0 best by absolute tau
        idx0 = sub["tau_level_0"].abs().idxmax()
        row0 = sub.loc[idx0]
        # Level 1 best by absolute tau
        idx1 = sub["tau_level_1"].abs().idxmax()
        row1 = sub.loc[idx1]
        rows.append({
            "file": param,
            "best_desc_level_0": row0["descriptor"],
            "tau_level_0": row0["tau_level_0"],
            "best_desc_level_1": row1["descriptor"],
            "tau_level_1": row1["tau_level_1"],
        })
    return pd.DataFrame(rows).sort_values("file")



def make_heatmap(matrix: pd.DataFrame, title: str, outpath: str):

    if HAS_SEABORN:
        sns.set_theme(style="white", context="paper", font_scale=1.2)

    fig, ax = plt.subplots(figsize=(14, 10))
    data = matrix.copy()

    data = data.reindex(sorted(data.columns), axis=1)
    vmin, vmax = -1.0, 1.0

    if HAS_SEABORN:
        from matplotlib.colors import BoundaryNorm
        import matplotlib.cm as cm

        boundaries = np.linspace(-1.0, 1.0, 9)
        norm = BoundaryNorm(boundaries, ncolors=256)

        sns.heatmap(
            data,
            cmap=CMAP,
            vmin=vmin,
            vmax=vmax,
            center=0,
            annot=False,
            cbar_kws={
                'label': 'Kendall τ Correlation',
                'shrink': 0.8,
                'aspect': 20,
                'pad': 0.02,
                'ticks': boundaries
            },
            linewidths=0.5,
            linecolor='white',
            square=False,
            cbar=True,
            ax=ax,
            norm=norm
        )

        ax.set_xlabel('flow param', fontsize=14, weight='bold')
        ax.set_ylabel('topo desc', fontsize=14, weight='bold')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    else:
        from matplotlib.colors import BoundaryNorm

        boundaries = np.linspace(-1.0, 1.0, 9)
        norm = BoundaryNorm(boundaries, ncolors=256)

        im = ax.imshow(data.values, aspect='auto', cmap=CMAP, vmin=vmin, vmax=vmax, norm=norm)
        cbar = plt.colorbar(im, ax=ax, label='Kendall τ Correlation', ticks=boundaries)
        ax.set_xticks(range(len(data.columns)))
        ax.set_xticklabels(data.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(data.index)))
        ax.set_yticklabels(data.index)
        ax.set_xlabel('flow param', fontsize=12, weight='bold')
        ax.set_ylabel('topo desc', fontsize=12, weight='bold')
    
    ax.set_title(title, fontsize=16, weight='bold', pad=20)
    
    plt.tight_layout()
    
    plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    if HAS_SEABORN:
        sns.reset_defaults()



os.makedirs(PAPER_RESULTS_DIR, exist_ok=True)
df = load_all_reports(RESULTS_DIR)

exclude_descriptors = ['FracD_x', 'FracD_y', 'FracD_k']
df_filtered = df[~df['descriptor'].isin(exclude_descriptors)].copy()

rank_df = compute_descriptor_rankings(df_filtered)
best_df = compute_best_descriptor_per_param(df_filtered)

mat0 = df_filtered.pivot_table(index="descriptor", columns="file", values="tau_level_0", aggfunc="first")
mat1 = df_filtered.pivot_table(index="descriptor", columns="file", values="tau_level_1", aggfunc="first")
order = rank_df["descriptor"].tolist()
mat0 = mat0.reindex(order)
mat1 = mat1.reindex(order)

make_heatmap(mat0, "tau heatmap L0", os.path.join(PAPER_RESULTS_DIR, "viz_heatmap_level_0.png"))
make_heatmap(mat1, "tau heatmap L1", os.path.join(PAPER_RESULTS_DIR, "viz_heatmap_level_1.png"))

rank_df.to_csv(os.path.join(PAPER_RESULTS_DIR, "viz_descriptor_rankings.csv"), index=False)
best_df.to_csv(os.path.join(PAPER_RESULTS_DIR, "viz_best_descriptors_by_parameter.csv"), index=False)

print("viz done")
