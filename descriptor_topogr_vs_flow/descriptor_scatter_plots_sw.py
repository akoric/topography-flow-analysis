import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def setup_publication_style():
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['axes.labelsize'] = 24
    plt.rcParams['axes.titlesize'] = 24
    plt.rcParams['xtick.labelsize'] = 22
    plt.rcParams['ytick.labelsize'] = 22
    plt.rcParams['legend.fontsize'] = 22
    plt.rcParams['font.size'] = 22
    plt.rcParams['lines.markersize'] = 8
    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Times']
    
    plt.rcParams['text.usetex'] = False  # Set to True if you have LaTeX installed
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    
    plt.rcParams['axes.grid'] = False
    plt.rcParams['axes.axisbelow'] = True
    
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.facecolor'] = 'white'

setup_publication_style()

def calculate_linear_regression(x, y):
    if len(x) < 2 or len(y) < 2:
        return None, None, 0.0
    
    x = np.array(x)
    y = np.array(y)
    
    coeffs = np.polyfit(x, y, 1)
    slope, intercept = coeffs
    
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
    
    return slope, intercept, r_squared

def load_roughness_data():
    roughness_data = {}
    
    for flow_param in ['U_mean', 'V_mean', 'W_mean']:
        file_path = f'data/roughness_{flow_param}.csv'
        try:
            df = pd.read_csv(file_path)
            # Extract DEM number from filename (e.g., DEM1_Moments.mat -> 1)
            df['DEM'] = df['File'].str.extract(r'DEM(\d+)').astype(int)
            roughness_data[flow_param] = df
        except FileNotFoundError:
            print(f"Warning: {file_path} not found")
            
    return roughness_data

def load_and_prepare_data_sw(topo_file='data/topography_results.csv', 
                            flow_file='data/U_mean_results.csv'):
    
    topo_data = pd.read_csv(topo_file)
    
    flow_data = pd.read_csv(flow_file)
    
    dem_frame_info = {}
    for dem in sorted(flow_data['DEM'].unique()):
        dem_frames = sorted(flow_data[flow_data['DEM'] == dem]['Frame'].unique())
        dem_frame_info[dem] = dem_frames
    
    return topo_data, flow_data, dem_frame_info

def create_multi_flow_scatter_plot_sw(descriptor, flow_params, save_path=None, figsize=(10, 8)):
    
    topo_file = f'data/topography_results.csv'
    topo_data = pd.read_csv(topo_file)
    
    if descriptor not in topo_data.columns:
        print(f"Descriptor '{descriptor}' not found in topography data")
        return None
    
    roughness_data = load_roughness_data()
    error_column = f'{descriptor}_std'
    
    topo_ordered = topo_data.sort_values(descriptor, ascending=True).reset_index(drop=True)
    rank_mapping = {int(row['DEM']): idx + 1 for idx, row in topo_ordered.iterrows()}
    
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    
    dem_markers = {
        1: 'o',      # circle
        2: 's',      # square
        3: '^',      # triangle up
        4: 'D',      # diamond
        5: 'v',      # triangle down
        6: 'p',      # pentagon
        7: '*',      # star
        8: 'h',      # hexagon
        9: 'X',      # X
        10: 'P'      # plus (filled)
    }
    
    all_data = {}
    
    for i, flow_param in enumerate(flow_params):
        flow_file = f'data/{flow_param}_results.csv'
        
        try:
            flow_data = pd.read_csv(flow_file)
        except FileNotFoundError:
            print(f"no file: {flow_file}")
            continue
            
        if descriptor not in flow_data.columns:
            print(f"Descriptor '{descriptor}' not found in {flow_param} data")
            continue
        
        dem_frame_info = {}
        for dem in sorted(flow_data['DEM'].unique()):
            dem_frames = sorted(flow_data[flow_data['DEM'] == dem]['Frame'].unique())
            dem_frame_info[dem] = dem_frames
        
        error_data = None
        if flow_param in roughness_data and error_column in roughness_data[flow_param].columns:
            error_data = roughness_data[flow_param]
        
        scatter_data = []
        
        for dem in topo_ordered['DEM']:
            dem_int = int(dem)
            topo_value = topo_data[topo_data['DEM'] == dem][descriptor].iloc[0]
            rank = rank_mapping[dem_int]
            
            # Get error bound for this DEM
            error_bound = 0
            if error_data is not None:
                error_row = error_data[error_data['DEM'] == dem_int]
                if len(error_row) > 0:
                    error_bound = error_row[error_column].iloc[0]
            
            if dem_int in dem_frame_info:
                frames = dem_frame_info[dem_int]
                if len(frames) >= 2:
                    level0_data = flow_data[(flow_data['DEM'] == dem_int) & 
                                          (flow_data['Frame'] == frames[0])]
                    if len(level0_data) > 0:
                        level0_value = level0_data[descriptor].iloc[0]
                        scatter_data.append({
                            'DEM': dem_int,
                            'rank': rank,
                            'topo_value': topo_value,
                            'flow_value': level0_value,
                            'error_bound': error_bound,
                            'level': 'Level 0',
                            'frame': frames[0],
                            'flow_param': flow_param
                        })
                    
                    level1_data = flow_data[(flow_data['DEM'] == dem_int) & 
                                          (flow_data['Frame'] == frames[1])]
                    if len(level1_data) > 0:
                        level1_value = level1_data[descriptor].iloc[0]
                        scatter_data.append({
                            'DEM': dem_int,
                            'rank': rank,
                            'topo_value': topo_value,
                            'flow_value': level1_value,
                            'error_bound': error_bound,
                            'level': 'Level 1', 
                            'frame': frames[1],
                            'flow_param': flow_param
                        })
        
        all_data[flow_param] = pd.DataFrame(scatter_data)
    
    if not all_data:
        print("No matching data found")
        return None
    
    fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
    
    for i, (flow_param, df) in enumerate(all_data.items()):
        level0_data = df[df['level'] == 'Level 0']
        if len(level0_data) > 0:
            # Plot each DEM separately to assign the correct marker
            for dem in sorted(level0_data['DEM'].unique()):
                dem_data = level0_data[level0_data['DEM'] == dem]
                marker = dem_markers.get(int(dem), 'o')
                
                if dem == sorted(level0_data['DEM'].unique())[0]:
                    label = f'{flow_param}'
                else:
                    label = None
                
                ax1.errorbar(dem_data['topo_value'], dem_data['flow_value'], 
                            yerr=dem_data['error_bound'],
                            fmt=marker, color=colors[i % len(colors)],
                            markersize=10, alpha=0.7, capsize=5, capthick=2,
                            label=label, linewidth=2, markeredgecolor='black', 
                            markeredgewidth=1.5)
    
    ax1.set_xlabel(f'{descriptor} topo')
    ax1.set_ylabel(f'{descriptor} flow')
    ax1.set_title(f'{descriptor} L0 comp')
    ax1.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        level0_path = save_path.replace('.png', '_level0.png')
        plt.savefig(level0_path, dpi=300, bbox_inches='tight')
    print(f"saved {level0_path}")
    
    plt.show()
    plt.close(fig1)
    
    fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
    
    for i, (flow_param, df) in enumerate(all_data.items()):
        level1_data = df[df['level'] == 'Level 1']
        if len(level1_data) > 0:
            # Plot each DEM separately to assign the correct marker
            for dem in sorted(level1_data['DEM'].unique()):
                dem_data = level1_data[level1_data['DEM'] == dem]
                marker = dem_markers.get(int(dem), 'o')
                
                if dem == sorted(level1_data['DEM'].unique())[0]:
                    label = f'{flow_param}'
                else:
                    label = None
                
                ax2.errorbar(dem_data['topo_value'], dem_data['flow_value'], 
                            yerr=dem_data['error_bound'],
                            fmt=marker, color=colors[i % len(colors)],
                            markersize=10, alpha=0.7, capsize=5, capthick=2,
                            label=label, linewidth=2, markeredgecolor='black', 
                            markeredgewidth=1.5)
    
    ax2.set_xlabel(f'{descriptor} topo')
    ax2.set_ylabel(f'{descriptor} flow')
    ax2.set_title(f'{descriptor} L1 comp')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        level1_path = save_path.replace('.png', '_level1.png')
        plt.savefig(level1_path, dpi=300, bbox_inches='tight')
    print(f"saved {level1_path}")
    
    plt.show()
    plt.close(fig2)
    
    return all_data

def create_level_comparison_plot_sw(descriptor, flow_params, save_path=None, figsize=(10, 8)):
    
    # Load topography data first
    topo_file = f'data/topography_results.csv'
    topo_data = pd.read_csv(topo_file)
    
    if descriptor not in topo_data.columns:
        print(f"Descriptor '{descriptor}' not found in topography data")
        return None
    
    # Load roughness data for error bounds
    roughness_data = load_roughness_data()
    error_column = f'{descriptor}_std'
    
    # Get topography ordering
    topo_ordered = topo_data.sort_values(descriptor, ascending=True).reset_index(drop=True)
    rank_mapping = {int(row['DEM']): idx + 1 for idx, row in topo_ordered.iterrows()}
    
    # Colors and markers for different flow parameters
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    markers = ['o', 's', '^', 'D', 'v']
    
    # Create the plot
    plt.figure(figsize=figsize)
    
    all_data = {}
    
    # Process each flow parameter
    for i, flow_param in enumerate(flow_params):
        flow_file = f'data/{flow_param}_results.csv'
        
        try:
            flow_data = pd.read_csv(flow_file)
        except FileNotFoundError:
            print(f"File not found: {flow_file}")
            continue
            
        if descriptor not in flow_data.columns:
            print(f"Descriptor '{descriptor}' not found in {flow_param} data")
            continue
        
        # Get frame info for this flow parameter
        dem_frame_info = {}
        for dem in sorted(flow_data['DEM'].unique()):
            dem_frames = sorted(flow_data[flow_data['DEM'] == dem]['Frame'].unique())
            dem_frame_info[dem] = dem_frames
        
        # Get error data for this flow parameter
        error_data = None
        if flow_param in roughness_data and error_column in roughness_data[flow_param].columns:
            error_data = roughness_data[flow_param]
        
        # Collect data for both levels
        scatter_data = []
        
        for dem in topo_ordered['DEM']:
            dem_int = int(dem)
            topo_value = topo_data[topo_data['DEM'] == dem][descriptor].iloc[0]
            rank = rank_mapping[dem_int]
            
            # Get error bound for this DEM
            error_bound = 0
            if error_data is not None:
                error_row = error_data[error_data['DEM'] == dem_int]
                if len(error_row) > 0:
                    error_bound = error_row[error_column].iloc[0]
            
            if dem_int in dem_frame_info:
                frames = dem_frame_info[dem_int]
                if len(frames) >= 2:
                    # Level 0: First frame
                    level0_data = flow_data[(flow_data['DEM'] == dem_int) & 
                                          (flow_data['Frame'] == frames[0])]
                    if len(level0_data) > 0:
                        level0_value = level0_data[descriptor].iloc[0]
                        scatter_data.append({
                            'DEM': dem_int,
                            'rank': rank,
                            'topo_value': topo_value,
                            'flow_value': level0_value,
                            'error_bound': error_bound,
                            'level': 'Level 0',
                            'frame': frames[0],
                            'flow_param': flow_param
                        })
                    
                    # Level 1: Second frame  
                    level1_data = flow_data[(flow_data['DEM'] == dem_int) & 
                                          (flow_data['Frame'] == frames[1])]
                    if len(level1_data) > 0:
                        level1_value = level1_data[descriptor].iloc[0]
                        scatter_data.append({
                            'DEM': dem_int,
                            'rank': rank,
                            'topo_value': topo_value,
                            'flow_value': level1_value,
                            'error_bound': error_bound,
                            'level': 'Level 1', 
                            'frame': frames[1],
                            'flow_param': flow_param
                        })
        
        all_data[flow_param] = pd.DataFrame(scatter_data)
    
    if not all_data:
        print("No matching data found")
        return None
    
    # Plot Level 0 vs Level 1 comparison with error bounds
    for i, (flow_param, df) in enumerate(all_data.items()):
        level0_data = df[df['level'] == 'Level 0']
        level1_data = df[df['level'] == 'Level 1']
        
        if len(level0_data) > 0 and len(level1_data) > 0:
            # Match DEMs between levels
            for dem in level0_data['DEM'].unique():
                l0_row = level0_data[level0_data['DEM'] == dem].iloc[0]
                l1_row = level1_data[level1_data['DEM'] == dem].iloc[0]
                
                l0_val = l0_row['flow_value']
                l1_val = l1_row['flow_value']
                error_bound = l0_row['error_bound']  # Same error for both levels
                
                # Plot point with error bars
                plt.errorbar(l0_val, l1_val, 
                           xerr=error_bound, yerr=error_bound,
                           fmt=markers[i % len(markers)], color=colors[i % len(colors)],
                           markersize=8, alpha=0.7, capsize=3, capthick=1,
                           label=f'{flow_param}' if dem == level0_data['DEM'].unique()[0] else "")
                
                # Use DEM number as label
                plt.annotate(f'{dem}', (l0_val, l1_val),
                            xytext=(3, 3), textcoords='offset points',
                            fontsize=10, alpha=0.9, fontweight='bold')
    
    # Add diagonal line for reference
    all_values = []
    for df in all_data.values():
        all_values.extend(df['flow_value'].tolist())
    if all_values:
        min_val, max_val = min(all_values), max(all_values)
        plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Equal values')
    
    plt.xlabel(f'{descriptor} - Level 0')
    plt.ylabel(f'{descriptor} - Level 1')
    plt.title(f'{descriptor}: Level 0 vs 1\n({", ".join(flow_params)})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Level comparison plot saved to: {save_path}")
    
    plt.show()
    
    return all_data

def create_ra_rq_sw_plots():
    
    Path('paper_results').mkdir(exist_ok=True)
    
    # Flow parameters to compare
    flow_params = ['U_mean', 'V_mean', 'W_mean']
    
    print(" Creating Ra_sw multi-flow scatter plot with error bounds")
    print("="*60)
    create_multi_flow_scatter_plot_sw('Ra_sw', flow_params, 
                                    save_path='paper_results/scatter_Ra_sw_U_mean_V_mean_W_mean.png')
    
    print("\n Creating Rq_sw multi-flow scatter plot with error bounds") 
    print("="*60)
    create_multi_flow_scatter_plot_sw('Rq_sw', flow_params,
                                    save_path='paper_results/scatter_Rq_sw_U_mean_V_mean_W_mean.png')
    
    print("\n Creating Ra_sw Level 0 vs Level 1 comparison plot with error bounds")
    print("="*70)
    create_level_comparison_plot_sw('Ra_sw', flow_params,
                                  save_path='paper_results/scatter_Ra_sw_level_comparison.png')
    
    print("\n Creating Rq_sw Level 0 vs Level 1 comparison plot with error bounds")
    print("="*70)
    create_level_comparison_plot_sw('Rq_sw', flow_params,
                                  save_path='paper_results/scatter_Rq_sw_level_comparison.png')

def demo_multi_flow_sw():
    create_ra_rq_sw_plots()

def quick_single_plot_sw(descriptor, flow_params=None):
    if flow_params is None:
        flow_params = ['U_mean', 'V_mean', 'W_mean']
    
    print(f" Creating {descriptor} multi-flow scatter plot with error bounds")
    print("="*60)
    save_path = f'paper_results/scatter_{descriptor}_{"_".join(flow_params)}.png'
    return create_multi_flow_scatter_plot_sw(descriptor, flow_params, save_path)


# Run demo
demo_multi_flow_sw()