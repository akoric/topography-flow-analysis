import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

DESCRIPTORS = [
    'FracD_x', 'Hurst_x', 'FracD_y', 'Hurst_y', 'FracD_k', 'Hurst_k',
    'Ra', 'Rq', 'Rq_slope', 'Sk', 'Ku', 'Lk', 'Lx', 'Ly',
    'Ra_sw', 'Rq_sw', 'Sk_sw', 'Ku_sw'
]


def find_csv_files(directory):
    pattern = os.path.join(directory, "*_results.csv")
    csv_files = glob.glob(pattern)
    return sorted([os.path.basename(f) for f in csv_files])

def load_and_align_data(csv_file_path):
    df = pd.read_csv(csv_file_path)
    
    dem_data = {}
    for dem in sorted(df['DEM'].unique()):
        dem_df = df[df['DEM'] == dem].copy()
        min_frame = dem_df['Frame'].min()
        dem_df['Aligned_Frame'] = dem_df['Frame'] - min_frame
        dem_df = dem_df.sort_values('Aligned_Frame')
        dem_data[int(dem)] = dem_df
    
    return dem_data


def create_descriptor_plots(dem_data, output_prefix, csv_filename):
    output_dir = f"plots_{output_prefix}"
    os.makedirs(output_dir, exist_ok=True)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for descriptor in DESCRIPTORS:
        plt.figure(figsize=(12, 8))
        
        for i, dem in enumerate(sorted(dem_data.keys())):
            df = dem_data[dem]
            if descriptor in df.columns:
                plt.plot(df['Aligned_Frame'], df[descriptor], 
                        color=colors[i], linewidth=2, marker='o', markersize=4,
                        label=f'DEM {dem}', alpha=0.8)
        
        plt.xlabel('Aligned Frame (starting from 0)', fontsize=12)
        plt.ylabel(descriptor, fontsize=12)
        plt.title(f'{descriptor} vs Frame\n({csv_filename})', fontsize=14, fontweight='bold')
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plot_filename = f"{output_dir}/{descriptor}_{output_prefix}.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()


def print_data_summary(dem_data, csv_filename):
    print(f"\nData Summary for {csv_filename}:")
    print("=" * 60)
    
    for dem in sorted(dem_data.keys()):
        df = dem_data[dem]
        original_min = df['Frame'].min()
        original_max = df['Frame'].max()
        aligned_max = df['Aligned_Frame'].max()
        n_frames = len(df)
        
        print(f"DEM {dem}:")
        print(f"  Original frames: {original_min:.0f} to {original_max:.0f}")
        print(f"  Aligned frames:  0 to {aligned_max:.0f}")
        print(f"  Total frames:    {n_frames}")
        print()


current_dir = "data"

csv_files = find_csv_files(current_dir)

if not csv_files:
    print("No _results.csv files found in the data directory!")
    exit()

print("Available CSV files:")
print("=" * 50)
for i, filename in enumerate(csv_files, 1):
    print(f"{i}. {filename}")

while True:
    try:
        choice = input(f"\nSelect a file to analyze (1-{len(csv_files)}) or 'q' to quit: ").strip()
        
        if choice.lower() == 'q':
            print("Goodbye!")
            exit()
        
        choice_idx = int(choice) - 1
        if 0 <= choice_idx < len(csv_files):
            selected_file = csv_files[choice_idx]
            break
        else:
            print(f"Please enter a number between 1 and {len(csv_files)}")
    except ValueError:
        print("Please enter a valid number or 'q' to quit")

print(f"\nSelected file: {selected_file}")

csv_path = os.path.join(current_dir, selected_file)

try:
    dem_data = load_and_align_data(csv_path)
    print_data_summary(dem_data, selected_file)
    output_prefix = selected_file.replace('_results.csv', '')
    create_descriptor_plots(dem_data, output_prefix, selected_file)
    print(f"\nDone. Plots in plots_{output_prefix}/")
    
except Exception as e:
    print(f"Error processing file: {e}")
    exit()
