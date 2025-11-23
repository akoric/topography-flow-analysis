import pandas as pd
import numpy as np
import os
import glob
from scipy.stats import kendalltau

def find_column_case_insensitive(df, possible_names):
    df_columns_lower = [col.lower() for col in df.columns]
    
    for name in possible_names:
        name_lower = name.lower()
        if name_lower in df_columns_lower:
            # Find the original column name
            idx = df_columns_lower.index(name_lower)
            return df.columns[idx]
    
    return None

def get_topography_order(topography_file, descriptor):
    print(f"\nuse topo file: {topography_file}")
    
    topo_df = pd.read_csv(topography_file)
    
    dem_col = find_column_case_insensitive(topo_df, ['DEM', 'dem', 'sample', 'id'])
    if dem_col is None:
        raise ValueError(f"Could not find DEM column in {topography_file}")
    
    if descriptor not in topo_df.columns:
        available_descriptors = [col for col in topo_df.columns if col not in [dem_col, 'Frame']]
        raise ValueError(f"Descriptor '{descriptor}' not found in topography file.\n"
                        f"Available descriptors: {available_descriptors}")
    
    dem_descriptor_pairs = []
    for _, row in topo_df.iterrows():
        dem = row[dem_col]
        desc_value = row[descriptor]
        if pd.notna(desc_value):  # Skip NaN values
            dem_descriptor_pairs.append((dem, desc_value))
    
    dem_descriptor_pairs.sort(key=lambda x: x[1])
    
    dem_order = {}
    for rank, (dem, _) in enumerate(dem_descriptor_pairs):
        dem_order[dem] = rank
    
    print(f"topo rank for '{descriptor}':")
    sorted_dems = sorted(dem_order.items(), key=lambda x: x[1])
    for dem, rank in sorted_dems:
        desc_val = topo_df[topo_df[dem_col] == dem][descriptor].iloc[0]
    print(f"  {rank}: DEM {dem} -> {desc_val:.6f}")
    
    return dem_order

def get_fluid_flow_orders(flow_file, descriptor, dem_order_reference):
    print(f"\nuse flow file: {flow_file}")
    
    flow_df = pd.read_csv(flow_file)
    
    dem_col = find_column_case_insensitive(flow_df, ['DEM', 'dem', 'sample', 'id'])
    frame_col = find_column_case_insensitive(flow_df, ['frame', 'z', 'level', 'Frame'])
    
    if dem_col is None:
        raise ValueError(f"Could not find DEM column in {flow_file}")
    if frame_col is None:
        raise ValueError(f"Could not find frame column in {flow_file}")
    
        if descriptor not in flow_df.columns:
            print(f"  warn: '{descriptor}' not in {flow_file}, skip")
            return None, None, 0, 0
    
    level_0_data = []
    level_1_data = []
    
    dem_groups = flow_df.groupby(dem_col)
    
    for dem, group in dem_groups:
        if dem not in dem_order_reference:
            continue
            
        sorted_frames = sorted(group[frame_col].unique())
        
        if len(sorted_frames) >= 1:
            level_0_frame = sorted_frames[0]
            level_0_rows = group[group[frame_col] == level_0_frame]
            if len(level_0_rows) > 0:
                desc_val = level_0_rows[descriptor].iloc[0]
                if pd.notna(desc_val):
                    level_0_data.append((dem, desc_val))
        
        if len(sorted_frames) >= 2:
            level_1_frame = sorted_frames[1]
            level_1_rows = group[group[frame_col] == level_1_frame]
            if len(level_1_rows) > 0:
                desc_val = level_1_rows[descriptor].iloc[0]
                if pd.notna(desc_val):
                    level_1_data.append((dem, desc_val))
    
    def create_ranking(data):
        if not data:
            return {}
        data.sort(key=lambda x: x[1])
        ranking = {}
        for rank, (dem, _) in enumerate(data):
            ranking[dem] = rank
        return ranking
    
    level_0_order = create_ranking(level_0_data)
    level_1_order = create_ranking(level_1_data)
    
    print(f"  L0 DEMs: {len(level_0_order)}")
    print(f"  L1 DEMs: {len(level_1_order)}")
    
    return level_0_order, level_1_order, len(level_0_order), len(level_1_order)

def compute_kendall_tau(order1, order2):
    common_dems = set(order1.keys()) & set(order2.keys())
    
    if len(common_dems) < 2:
        return np.nan
    
    ranks1 = [order1[dem] for dem in common_dems]
    ranks2 = [order2[dem] for dem in common_dems]
    
    tau, p_value = kendalltau(ranks1, ranks2)
    
    return tau

data_dir = "data"
results_dir = "results"
topography_file = os.path.join(data_dir, "topography_results.csv")
output_file = os.path.join(results_dir, "kendall_report.csv")

os.makedirs(results_dir, exist_ok=True)

if not os.path.exists(data_dir):
    print(f"err: data dir '{data_dir}' not found")
    exit(1)

if not os.path.exists(topography_file):
    print(f"err: topo file '{topography_file}' not found")
    exit(1)

topo_df = pd.read_csv(topography_file)
dem_col = find_column_case_insensitive(topo_df, ['DEM', 'dem', 'sample', 'id'])
available_descriptors = [col for col in topo_df.columns if col not in [dem_col, 'Frame']]

print("kendall correlation tool")
print(f"descriptors: {available_descriptors}")

while True:
    descriptor = input(f"desc name: ").strip()
    if descriptor in available_descriptors:
        break
    print("bad desc, try again")

print(f"\nuse desc: {descriptor}")

try:
    topo_order = get_topography_order(topography_file, descriptor)

    flow_files = glob.glob(os.path.join(data_dir, "*.csv"))
    flow_files = [f for f in flow_files if not f.endswith("topography_results.csv")]
    flow_files = [f for f in flow_files if not os.path.basename(f).startswith("kendall_report")]
    flow_files.sort()

    print(f"\nfound {len(flow_files)} flow files")

    results = []

    for flow_file in flow_files:
        file_name = os.path.basename(flow_file)
        print(f"\nfile: {file_name}")

        try:
            level_0_order, level_1_order, n_0, n_1 = get_fluid_flow_orders(
                flow_file, descriptor, topo_order)

            if level_0_order is None:
                continue

            tau_0 = compute_kendall_tau(topo_order, level_0_order)
            tau_1 = compute_kendall_tau(topo_order, level_1_order)

            results.append({
                'file': file_name,
                'tau_level_0': tau_0,
                'tau_level_1': tau_1,
                'n_DEM_level_0': n_0,
                'n_DEM_level_1': n_1
            })

            if not np.isnan(tau_0):
                print(f"  tau L0: {tau_0:.4f}")
            else:
                print("  tau L0: N/A")

            if not np.isnan(tau_1):
                print(f"  tau L1: {tau_1:.4f}")
            else:
                print("  tau L1: N/A")

        except Exception as e:
            print(f"  err {file_name}: {e}")
            continue

    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, index=False)

        print("\nsummary:")
        print(results_df.to_string(index=False, float_format='%.4f'))
        print(f"saved {output_file}")
    else:
        print("\nno results (errors or missing desc)")

except Exception as e:
    print(f"\nfatal err: {e}")
    exit(1)
