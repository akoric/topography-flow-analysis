import os
import pandas as pd
import subprocess
import sys

def get_available_descriptors():
    topography_file = "data/topography_results.csv"
    
    if not os.path.exists(topography_file):
        print(f"Error: {topography_file} not found")
        return []
    
    df = pd.read_csv(topography_file)
    
    dem_cols = ['DEM', 'dem', 'sample', 'id']
    dem_col = None
    for col in dem_cols:
        if col in df.columns:
            dem_col = col
            break
    
    if dem_col is None:
        print("Error: Could not find DEM column")
        return []
    
    descriptors = [col for col in df.columns if col not in [dem_col, 'Frame']]
    return descriptors

def run_analysis_for_descriptor(descriptor):
    print(f"\nrun {descriptor}")

    try:
        result = subprocess.run(
            ['python3', 'kendall_correlation_analysis.py'],
            input=descriptor + '\n',
            text=True,
            capture_output=True
        )
        
        if result.returncode == 0:
            print(result.stdout)

            original_file = "results/kendall_report.csv"
            new_file = f"results/kendall_report_{descriptor}.csv"
            
            if os.path.exists(original_file):
                os.rename(original_file, new_file)
            
            return True
        else:
            print(f"Error running analysis for {descriptor}:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"Error running {descriptor}: {e}")
        return False
print("batch kendall analysis")

descriptors = get_available_descriptors()

if not descriptors:
    print("no descriptors, exit")
    sys.exit(1)

print(f"descriptors: {descriptors}")

print("1: all, 2: pick, 3: list")
choice = input("choice: ").strip()

if choice == "3":
    print("available:")
    for i, desc in enumerate(descriptors, 1):
        print(f"{i:2d}: {desc}")
    sys.exit(0)
elif choice == "1":
    selected_descriptors = descriptors
elif choice == "2":
    print("available:")
    for i, desc in enumerate(descriptors, 1):
        print(f"{i:2d}: {desc}")

    selection = input("nums (comma): ").strip()

    try:
        indices = [int(x.strip()) - 1 for x in selection.split(',')]
        selected_descriptors = [descriptors[i] for i in indices if 0 <= i < len(descriptors)]

        if not selected_descriptors:
            print("no valid selection")
            sys.exit(1)

    except (ValueError, IndexError):
        print("bad selection")
        sys.exit(1)
else:
    print("bad choice")
    sys.exit(1)

print(f"run {len(selected_descriptors)}: {selected_descriptors}")

try:
    confirm = input("go? (y/n): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("cancelled")
        sys.exit(0)
except EOFError:
    confirm = 'y'

successful = 0
failed = 0

for descriptor in selected_descriptors:
    if run_analysis_for_descriptor(descriptor):
        successful += 1
    else:
        failed += 1

print("done batch")
print(f"ok: {successful}")
print(f"fail: {failed}")
print(f"total: {len(selected_descriptors)}")
