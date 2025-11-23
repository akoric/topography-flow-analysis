import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os

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

def get_available_descriptors():
    topo_file = 'data/topography_results.csv'
    if not os.path.exists(topo_file):
        print(f"Error: {topo_file} not found")
        return []
    
    topo_data = pd.read_csv(topo_file)
    exclude_cols = ['DEM', 'Frame', 'dem', 'sample', 'id']
    descriptors = [col for col in topo_data.columns if col not in exclude_cols]
    return sorted(descriptors)

def get_available_flow_parameters():
    data_dir = 'data'
    if not os.path.exists(data_dir):
        print(f"Error: {data_dir} not found")
        return []
    
    # Find all CSV files matching pattern: {direction}_{param}_results.csv
    import glob
    files = glob.glob(os.path.join(data_dir, '*_results.csv'))
    
    # Extract unique parameter types
    param_types = set()
    for file in files:
        basename = os.path.basename(file)
        # Remove '_results.csv'
        name = basename.replace('_results.csv', '')
        # Check if it starts with U_, V_, or W_
        if name.startswith('U_') or name.startswith('V_') or name.startswith('W_'):
            param = name[2:]  # Remove 'U_', 'V_', 'W_' prefix
            param_types.add(param)
    
    return sorted(list(param_types))

def create_multi_direction_scatter_plot(descriptor, flow_param_type, save_path=None, figsize=(10, 8)):

    subdir_name = f'{descriptor}_vs_{flow_param_type}'
    output_dir = Path('paper_results') / subdir_name
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"dir: {output_dir}")
    
    # Load topography data first
    topo_file = f'data/topography_results.csv'
    topo_data = pd.read_csv(topo_file)
    
    if descriptor not in topo_data.columns:
        print(f"Descriptor '{descriptor}' not found in topography data")
        return None
    
    # Get topography ordering
    topo_ordered = topo_data.sort_values(descriptor, ascending=True).reset_index(drop=True)
    rank_mapping = {int(row['DEM']): idx + 1 for idx, row in topo_ordered.iterrows()}
    
    colors = ['blue', 'green', 'red']
    directions = ['U', 'V', 'W']
    
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
    
    for i, direction in enumerate(directions):
        flow_param = f'{direction}_{flow_param_type}'
        flow_file = f'data/{flow_param}_results.csv'
        
        if not os.path.exists(flow_file):
            print(f"warn: no {flow_file}")
            continue
            
        flow_data = pd.read_csv(flow_file)
        
        if descriptor not in flow_data.columns:
            print(f"warn: {descriptor} not in {flow_param}")
            continue
        
        dem_frame_info = {}
        for dem in sorted(flow_data['DEM'].unique()):
            dem_frames = sorted(flow_data[flow_data['DEM'] == dem]['Frame'].unique())
            dem_frame_info[dem] = dem_frames
        
        scatter_data = []
        
        for dem in topo_ordered['DEM']:
            dem_int = int(dem)
            topo_value = topo_data[topo_data['DEM'] == dem][descriptor].iloc[0]
            rank = rank_mapping[dem_int]
            
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
                            'level': 'Level 0',
                            'frame': frames[0],
                            'direction': direction
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
                            'level': 'Level 1', 
                            'frame': frames[1],
                            'direction': direction
                        })
        
        all_data[direction] = pd.DataFrame(scatter_data)
    
    if not all_data:
        print("No matching data found")
        return None
    
    fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
    
    for i, (direction, df) in enumerate(all_data.items()):
        level0_data = df[df['level'] == 'Level 0']
        if len(level0_data) > 0:
            # Plot each DEM separately to assign the correct marker
            for dem in sorted(level0_data['DEM'].unique()):
                dem_data = level0_data[level0_data['DEM'] == dem]
                marker = dem_markers.get(int(dem), 'o')  # Default to circle if DEM not in dict
                
                if dem == sorted(level0_data['DEM'].unique())[0]:
                    label = f'{direction}_{flow_param_type}'
                else:
                    label = None
                
                ax1.scatter(dem_data['topo_value'], dem_data['flow_value'], 
                           c=colors[i % len(colors)], s=150, alpha=0.7, 
                           label=label, marker=marker, edgecolors='black', linewidths=1.5)
    
    ax1.set_xlabel(f'{descriptor} topo')
    ax1.set_ylabel(f'{descriptor} flow')
    ax1.set_title(f'{descriptor} vs {flow_param_type.upper()} L0 (U,V,W)')
    ax1.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        level0_path = output_dir / f'scatter_{descriptor}_{flow_param_type}_UVW_level0.png'
        plt.savefig(level0_path, dpi=300, bbox_inches='tight')
    print(f"saved {level0_path}")
    
    plt.show()
    plt.close(fig1)
    
    fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
    
    for i, (direction, df) in enumerate(all_data.items()):
        level1_data = df[df['level'] == 'Level 1']
        if len(level1_data) > 0:
            # Plot each DEM separately to assign the correct marker
            for dem in sorted(level1_data['DEM'].unique()):
                dem_data = level1_data[level1_data['DEM'] == dem]
                marker = dem_markers.get(int(dem), 'o')  # Default to circle if DEM not in dict
                
                if dem == sorted(level1_data['DEM'].unique())[0]:
                    label = f'{direction}_{flow_param_type}'
                else:
                    label = None
                
                ax2.scatter(dem_data['topo_value'], dem_data['flow_value'], 
                           c=colors[i % len(colors)], s=150, alpha=0.7, 
                           label=label, marker=marker, edgecolors='black', linewidths=1.5)
    
    ax2.set_xlabel(f'{descriptor} topo')
    ax2.set_ylabel(f'{descriptor} flow')
    ax2.set_title(f'{descriptor} vs {flow_param_type.upper()} L1 (U,V,W)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        level1_path = output_dir / f'scatter_{descriptor}_{flow_param_type}_UVW_level1.png'
        plt.savefig(level1_path, dpi=300, bbox_inches='tight')
    print(f"saved {level1_path}")
    
    plt.show()
    plt.close(fig2)
    
    return all_data

def interactive_mode():
    print("scatter plot tool")
    
    # Get available options
    descriptors = get_available_descriptors()
    flow_params = get_available_flow_parameters()
    
    if not descriptors:
        print("no topo descriptors")
        return
    
    if not flow_params:
        print("no flow params")
        return
    
    print("\ntopo desc:")
    for i, desc in enumerate(descriptors, 1):
        print(f"  {i:2d}. {desc}")
    
    # Get descriptor choice
    while True:
        try:
            choice = input(f"\nEnter descriptor number (1-{len(descriptors)}) or name: ").strip()
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(descriptors):
                    descriptor = descriptors[idx]
                    break
            elif choice in descriptors:
                descriptor = choice
                break
            print("bad choice")
        except (ValueError, KeyboardInterrupt):
            print("\nOperation cancelled.")
            return
    
    print(f"\nuse desc: {descriptor}")
    
    print("\nflow types:")
    for i, param in enumerate(flow_params, 1):
        print(f"  {i:2d}. {param}")
    
    # Get flow parameter choice
    while True:
        try:
            choice = input(f"\nEnter flow parameter number (1-{len(flow_params)}) or name: ").strip()
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(flow_params):
                    flow_param = flow_params[idx]
                    break
            elif choice in flow_params:
                flow_param = choice
                break
            print("bad choice")
        except (ValueError, KeyboardInterrupt):
            print("\nOperation cancelled.")
            return
    
    print(f"\nuse flow: {flow_param}")

    create_multi_direction_scatter_plot(descriptor, flow_param, save_path='dummy')
    print("done")

def quick_plot(descriptor, flow_param_type):
    return create_multi_direction_scatter_plot(descriptor, flow_param_type, save_path='dummy')

interactive_mode()
