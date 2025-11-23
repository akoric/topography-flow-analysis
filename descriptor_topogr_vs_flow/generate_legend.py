import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Times']
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['savefig.facecolor'] = 'white'

def create_compact_legend():

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    dem_markers = {1: 'o', 2: 's', 3: '^', 4: 'D'}
    
    ax1.axis('off')
    ax1.set_title('DEM', fontsize=22, fontweight='bold')
    
    dem_legend = []
    for dem in sorted(dem_markers.keys()):
        dem_legend.append(
            Line2D([0], [0], marker=dem_markers[dem], color='w', 
                   label=f'DEM {dem}',
                   markerfacecolor='gray', markersize=18,
                   markeredgecolor='black', markeredgewidth=2)
        )
    
    ax1.legend(handles=dem_legend, loc='center', frameon=True,
              fancybox=True, shadow=True, fontsize=20)
    
    ax2.axis('off')
    ax2.set_title('Flow', fontsize=22, fontweight='bold')
    
    flow_colors = {'U': 'blue', 'V': 'green', 'W': 'red'}
    flow_legend = []
    for direction, color in flow_colors.items():
        flow_legend.append(
            Line2D([0], [0], marker='o', color='w',
                   label=direction,
                   markerfacecolor=color, markersize=18,
                   markeredgecolor='black', markeredgewidth=2)
        )
    
    ax2.legend(handles=flow_legend, loc='center', frameon=True,
              fancybox=True, shadow=True, fontsize=20)
    
    plt.suptitle('Legend', fontsize=24, fontweight='bold')
    plt.tight_layout()
    
    output_dir = Path('paper_results')
    output_dir.mkdir(exist_ok=True)
    save_path = output_dir / 'legend_compact.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"saved {save_path}")
    
    plt.show()
    plt.close()

create_compact_legend()