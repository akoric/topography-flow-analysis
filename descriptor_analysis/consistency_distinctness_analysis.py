import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os
import glob
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
sns.set_palette("husl")

DESCRIPTORS = [
    'FracD_x', 'Hurst_x', 'FracD_y', 'Hurst_y', 'FracD_k', 'Hurst_k',
    'Ra', 'Rq', 'Rq_slope', 'Sk', 'Ku', 'Lk', 'Lx', 'Ly',
    'Ra_sw', 'Rq_sw', 'Sk_sw', 'Ku_sw'
]


class ImprovedSlopeInterceptAnalyzer:
    
    def __init__(self, num_frames=4):
        self.num_frames = num_frames
        self.normalize_method = 'z'
        self.normalize_scope = 'frame0'
        self.log1p = False
        self.results = {}
        
    def load_and_prepare_data(self, csv_file_path):
        df = pd.read_csv(csv_file_path)
        
        aligned_data = []
        for dem in sorted(df['DEM'].unique()):
            dem_df = df[df['DEM'] == dem].copy()
            min_frame = dem_df['Frame'].min()
            dem_df['Aligned_Frame'] = dem_df['Frame'] - min_frame
            aligned_data.append(dem_df)
        
        self.data = pd.concat(aligned_data, ignore_index=True)
        self.data = self.data[self.data['Aligned_Frame'] < self.num_frames]
        self._normalize_descriptors_in_place()
        
        return self.data
    
    def _normalize_descriptors_in_place(self):
        eps = 1e-12
        for desc in DESCRIPTORS:
            if desc not in self.data.columns:
                continue

            ref = self.data[self.data['Aligned_Frame'] == 0][desc].astype(float)

            if ref.size < 2:
                continue

            x = self.data[desc].astype(float).values

            mu = float(ref.mean())
            sd = float(ref.std(ddof=1))
            if sd <= 0:
                continue
            x_norm = (x - mu) / (sd + eps)

            self.data[desc] = x_norm
    
    def _interpret_directional_consistency(self, consistency_score, mean_slope):
        if consistency_score >= 0.8:
            consistency_level = "very high"
        elif consistency_score >= 0.6:
            consistency_level = "high"
        elif consistency_score >= 0.4:
            consistency_level = "moderate"
        elif consistency_score >= 0.2:
            consistency_level = "low"
        else:
            consistency_level = "very low"
        
        abs_mean_slope = abs(mean_slope)
        if abs_mean_slope < 0.01:
            trend = "flat"
        elif mean_slope > 0:
            trend = "upward" if abs_mean_slope <= 0.1 else "strongly upward"
        else:
            trend = "downward" if abs_mean_slope <= 0.1 else "strongly downward"
        
        if consistency_score >= 0.6:
            return f"Trend: {trend} ({consistency_level} consistency)"
        else:
            return f"Mixed trends ({consistency_level} consistency)"
    
    def analyze_descriptor(self, descriptor):
        if descriptor not in self.data.columns:
            return None
            
        desc_data = self.data[['DEM', 'Aligned_Frame', descriptor]].dropna()
        
        if len(desc_data) < 10:
            return None
            
        desc_data['DEM_cat'] = desc_data['DEM'].astype(str)
        
        results = {
            'descriptor': descriptor,
            'n_observations': len(desc_data),
            'dems': sorted(desc_data['DEM'].unique())
        }
        
        frame_separations = []
        all_frame_intercepts = {}
        
        for frame in range(self.num_frames):
            frame_data = desc_data[desc_data['Aligned_Frame'] == frame]
            if len(frame_data) >= len(results['dems']):
                try:
                    frame_values = []
                    for dem in results['dems']:
                        dem_frame = frame_data[frame_data['DEM'] == dem][descriptor].values
                        if len(dem_frame) > 0:
                            frame_values.append(float(dem_frame[0]))
                    
                    all_frame_intercepts[frame] = frame_values
                    
                    if len(frame_values) >= 2:
                        pairwise_distances = []
                        for i in range(len(frame_values)):
                            for j in range(i + 1, len(frame_values)):
                                distance = abs(frame_values[i] - frame_values[j])
                                pairwise_distances.append(distance)
                        
                        min_separation = min(pairwise_distances) if pairwise_distances else 0
                        frame_separations.append(min_separation)
                        
                except:
                    continue
        
        if frame_separations:
            results['min_separation_across_frames'] = min(frame_separations)
            results['mean_separation_across_frames'] = np.mean(frame_separations)
            results['frame_separations'] = frame_separations
        else:
            results['min_separation_across_frames'] = 0
            results['mean_separation_across_frames'] = 0
            results['frame_separations'] = []
        
        dem_fits = {}
        slopes = []
        intercepts = []
        r_squares = []
        
        for dem in results['dems']:
            dem_data = desc_data[desc_data['DEM'] == dem]
            if len(dem_data) >= 2:
                try:
                    X = dem_data['Aligned_Frame'].values
                    y = dem_data[descriptor].values
                    
                    X_with_const = sm.add_constant(X)
                    model = sm.OLS(y, X_with_const).fit()
                    
                    intercept = model.params[0]
                    slope = model.params[1]
                    
                    dem_fits[dem] = {
                        'intercept': intercept,
                        'slope': slope,
                        'r_squared': model.rsquared,
                        'model': model
                    }
                    
                    slopes.append(slope)
                    intercepts.append(intercept)
                    r_squares.append(model.rsquared)
                    
                except:
                    continue
        
        results['dem_fits'] = dem_fits
        results['slopes'] = slopes
        results['intercepts'] = intercepts
        
        if slopes:
            results['mean_slope'] = np.mean(slopes)
            results['slope_std'] = np.std(slopes, ddof=1) if len(slopes) > 1 else 0
            
            slope_std_dev = np.std(slopes, ddof=1) if len(slopes) > 1 else 0
            results['directional_consistency'] = 1.0 / (1.0 + slope_std_dev)
            results['slope_stability'] = results['directional_consistency']
        else:
            results['mean_slope'] = 0
            results['slope_std'] = 0
            results['directional_consistency'] = 0
            results['slope_stability'] = 0
        
        results['mean_r_squared'] = np.mean(r_squares) if r_squares else 0
        
        if results['mean_separation_across_frames'] is not None and results['slope_stability'] is not None:
            results['combined_score'] = results['mean_separation_across_frames'] * results['slope_stability']
        else:
            results['combined_score'] = 0.0
        
        results['meets_criteria'] = len(slopes) >= 2
        
        return results
    
    def analyze_all_descriptors(self):
        all_results = []
        
        for descriptor in DESCRIPTORS:
            result = self.analyze_descriptor(descriptor)
            if result is not None:
                all_results.append(result)
        
        all_results.sort(key=lambda x: (not x['meets_criteria'], -x['combined_score']))
        
        self.results = all_results
        return all_results
    
    def print_summary(self):
        if not self.results:
            print("No results to display. Run analyze_all_descriptors() first.")
            return
        
        print(f"\n{'='*90}")
        print("ANALYSIS SUMMARY")
        print(f"{'='*90}")
        
        meeting_criteria = [r for r in self.results if r['meets_criteria']]
        print(f"Valid descriptors: {len(meeting_criteria)}/{len(self.results)}")
        print()
        
        if meeting_criteria:
            print("TOP DESCRIPTORS:")
            print("-" * 90)
            print(f"{'Rank':<4} {'Descriptor':<12} {'Avg Sep':<8} {'Slope SD':<10} {'Consistency':<11} {'Combined':<10} {'R²':<6}")
            print("-" * 90)
            
            for i, result in enumerate(meeting_criteria, 1):
                print(f"{i:<4} {result['descriptor']:<12} {result.get('mean_separation_across_frames', 0):<8.3f} "
                      f"{result['slope_std']:<10.3f} {result['directional_consistency']:<11.3f} "
                      f"{result['combined_score']:<10.3f} {result['mean_r_squared']:<6.3f}")
        
        print("\nALL DESCRIPTORS:")
        print("-" * 80)
        print(f"{'Rank':<4} {'Descriptor':<12} {'Avg Sep':<8} {'Slope SD':<10} {'Combined':<10}")
        print("-" * 80)
        
        for i, result in enumerate(self.results, 1):
            print(f"{i:<4} {result['descriptor']:<12} {result.get('mean_separation_across_frames', 0):<8.3f} "
                  f"{result['slope_std']:<10.3f} {result['combined_score']:<10.3f}")
    
    def create_plots(self, output_dir="improved_slope_intercept_plots", max_plots=10):
        if not self.results:
            return
        
        os.makedirs(output_dir, exist_ok=True)
        descriptors_to_plot = self.results[:max_plots]
        
        for i, result in enumerate(descriptors_to_plot, 1):
            self._create_descriptor_plot(result, output_dir, rank=i)
    
    def _create_descriptor_plot(self, result, output_dir, rank=None):
        """Create a detailed plot for a single descriptor"""
        descriptor = result['descriptor']
        
        # prepare data
        desc_data = self.data[['DEM', 'Aligned_Frame', descriptor]].dropna()
        
        # create figure with subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot 1: Data points and fitted lines
        colors = plt.cm.Set1(np.linspace(0, 1, len(result['dems'])))
        
        for dem_idx, (dem, color) in enumerate(zip(result['dems'], colors)):
            dem_data = desc_data[desc_data['DEM'] == dem]
            
            # plot data points
            ax1.scatter(dem_data['Aligned_Frame'], dem_data[descriptor], 
                       color=color, alpha=0.7, s=50, label=f'DEM {dem}')
            
            # plot fitted line if available
            if dem in result['dem_fits']:
                fit = result['dem_fits'][dem]
                x_line = np.linspace(0, self.num_frames-1, 100)
                y_line = fit['intercept'] + fit['slope'] * x_line
                ax1.plot(x_line, y_line, color=color, linewidth=2, alpha=0.8,
                        linestyle='--')
                
                # add slope annotation using enumeration index to prevent overflow
                ax1.text(0.02, 0.98 - dem_idx*0.05, f'DEM {dem}: slope = {fit["slope"]:.4f}',
                        transform=ax1.transAxes, fontsize=9, color=color,
                        verticalalignment='top')
        
        ax1.set_xlabel('Aligned Frame')
        ax1.set_ylabel(descriptor)
        ax1.set_title(f'{descriptor} vs Frame (Rank #{rank if rank else "N/A"})')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Improved statistics summary
        ax2.axis('off')
        
        stats_text = f"""
ANALYSIS RESULTS
{'-'*35}

Separation:
• Avg: {result.get('mean_separation_across_frames', 0):.3f} SD
• Min: {result.get('min_separation_across_frames', 0):.3f} SD

Slope Stats:
• Std Dev: {result['slope_std']:.4f}
• Mean: {result['mean_slope']:.4f}
• Consistency: {result['directional_consistency']:.3f}

{self._interpret_directional_consistency(result['directional_consistency'], result['mean_slope'])}

Combined: {result['combined_score']:.3f}
R²: {result['mean_r_squared']:.3f}
N: {result['n_observations']}
        """
        
        ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes, 
                fontsize=7, verticalalignment='top', fontfamily='monospace')
        
        if result['frame_separations']:
            ax2_inset = fig.add_axes([0.55, 0.02, 0.35, 0.15])
            frames = list(range(len(result['frame_separations'])))
            ax2_inset.bar(frames, result['frame_separations'], alpha=0.7, color='skyblue')
            ax2_inset.set_xlabel('Frame', fontsize=7)
            ax2_inset.set_ylabel('Min Sep', fontsize=7)
            ax2_inset.set_title('Separation by Frame', fontsize=8)
            ax2_inset.tick_params(labelsize=6)
        
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.20, wspace=0.3)
        
        if rank is not None:
            plot_filename = f"{output_dir}/{rank:02d}_{descriptor}_improved_analysis.png"
        else:
            plot_filename = f"{output_dir}/{descriptor}_improved_analysis.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
    
    def save_results_csv(self, output_file="improved_slope_intercept_results.csv"):
        if not self.results:
            return
        
        csv_data = []
        for result in self.results:
            row = {
                'descriptor': result['descriptor'],
                'mean_separation_across_frames': result.get('mean_separation_across_frames', 0),
                'min_separation_across_frames': result.get('min_separation_across_frames', 0),
                'directional_consistency': result.get('directional_consistency', 0),
                'slope_stability': result['slope_stability'],
                'combined_score': result['combined_score'],
                'mean_slope': result['mean_slope'],
                'slope_std': result['slope_std'],
                'mean_r_squared': result['mean_r_squared'],
                'n_observations': result['n_observations']
            }
            
            for i, sep in enumerate(result.get('frame_separations', [])):
                row[f'frame_{i}_separation'] = sep
            
            for dem in result['dems']:
                if dem in result['dem_fits']:
                    fit = result['dem_fits'][dem]
                    row[f'DEM_{dem}_intercept'] = fit['intercept']
                    row[f'DEM_{dem}_slope'] = fit['slope']
                    row[f'DEM_{dem}_r_squared'] = fit['r_squared']
            
            csv_data.append(row)
        
        df_results = pd.DataFrame(csv_data)
        df_results.to_csv(output_file, index=False)


def find_csv_files(directory):
    pattern = os.path.join(directory, "*_results.csv")
    csv_files = glob.glob(pattern)
    return sorted([os.path.basename(f) for f in csv_files])


data_dir = "data"

csv_files = find_csv_files(data_dir)

if not csv_files:
    print("No _results.csv files found in the data directory!")
    exit()

print("Available CSV files:")
print("=" * 50)
for i, filename in enumerate(csv_files, 1):
    print(f"{i:2d}. {filename}")

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
csv_path = os.path.join(data_dir, selected_file)

try:
    num_frames = int(input("Number of frames (default 4): ").strip() or "4")
except ValueError:
    num_frames = 4

analyzer = ImprovedSlopeInterceptAnalyzer(num_frames=num_frames)

try:
    analyzer.load_and_prepare_data(csv_path)
    analyzer.analyze_all_descriptors()
    analyzer.print_summary()
    
    output_prefix = selected_file.replace('_results.csv', '')
    plot_dir = f"plots_{output_prefix}_consistency_distinctness_analysis"
    analyzer.create_plots(output_dir=plot_dir, max_plots=18)
    
    results_file = f"{output_prefix}_consistency_distinctness_analysis.csv"
    analyzer.save_results_csv(results_file)
    
    print(f"\nDone. Results in {results_file}")
    
except Exception as e:
    print(f"Error: {e}")
