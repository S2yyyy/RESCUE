"""Shared plotting style for the RESCUE paper."""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({
    'font.size': 9,
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'text.usetex': False,
    'mathtext.fontset': 'stix',
    'lines.linewidth': 1.2,
    'patch.linewidth': 0.6,
})

COLORS = {
    'blue': '#2166AC',
    'red': '#B2182B',
    'green': '#4DAF4A',
    'orange': '#FF7F00',
    'purple': '#984EA3',
    'gray': '#666666',
    'light_gray': '#BDBDBD',
}

SPECIES_COLORS = {
    'mouse': '#2166AC',
    'dog': '#4DAF4A',
    'chicken': '#FF7F00',
    'zebrafish': '#B2182B',
}

SPECIES_LABELS = {
    'mouse': 'Mouse',
    'dog': 'Dog',
    'chicken': 'Chicken',
    'zebrafish': 'Zebrafish',
}

TOOL_COLORS = {
    'miniprot': '#2166AC',
    'Spaln3': '#B2182B',
    'GTH': '#4DAF4A',
    'Scipio': '#FF7F00',
    'genBlastG': '#984EA3',
    'mini+miniSplice': '#A65628',
    'Genewise': '#F781BF',
    'Exonerate': '#999999',
}

import os
FIG_DIR = os.path.dirname(os.path.abspath(__file__))

def save_fig(fig, name):
    path = f'{FIG_DIR}/{name}.pdf'
    fig.savefig(path)
    print(f'Saved: {path}')
    plt.close(fig)
