import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from scipy import stats


allpass = pd.read_csv("sheafhasshetero.csv")

colors = {'b': "#50DD71",
          'n': "#65ADFF",
          'o': "#FF5050",  
          's': '#FF9500'}   

labels = {'b': 'Boron', 's': 'Sulfur', 'o': 'Oxygen', 'n': 'Nitrogen'}

# heteroatoms combined plot
allpass['color'] = (allpass[['b', 's', 'o', 'n']].idxmax(axis=1).map(colors))

plt.figure(figsize=(10, 6))
plt.xlabel('DFT HOMO-LUMO gap / eV')
plt.ylabel('Frontier eigengap')
plt.grid(visible=True, which="both", color="0.8", zorder=1)

for atom, color in colors.items():
    subset = allpass[allpass[atom] > 0]
    plt.scatter(subset["gap"], subset["sheafGAP"], facecolor=color, edgecolor=color, alpha=0.7, marker=".", s=2, zorder=2)

x = allpass['gap']
y = allpass['sheafGAP']

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
x_line = np.linspace(x.min(), x.max(), 100)
plt.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)

plt.text(0.05, 0.95, f'y = {slope:.3f}x {intercept:+.3f}\n$R^2$ = {r_value**2:.3f}', transform=plt.gca().transAxes, verticalalignment='top', fontsize=10)
legend = [mpatches.Patch(color=color, alpha=0.7, label=labels[atom]) for atom, color in colors.items()]
plt.legend(handles=legend, title='Heteroatom', markerscale=3)

plt.xlim(left=0)
plt.ylim(bottom=0)

plt.savefig("HASgaphetero.pdf", format="pdf", bbox_inches='tight')
plt.clf()


# seaprated by heteroatom


subsets = {'b': allpass[allpass['b'] > 0],
           's': allpass[allpass['s'] > 0],
           'o': allpass[allpass['o'] > 0],
           'n': allpass[allpass['n'] > 0]}

for atom, color in colors.items():
    label = labels[atom]
    subset = subsets[atom]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_xlabel('DFT HOMO-LUMO gap / eV')
    ax.set_ylabel('Frontier eigengap')
    ax.set_title(f'{label}-containing molecules (n={len(subset)})')
    ax.grid(visible=True, which="both", color="0.8", zorder=1)
    ax.scatter(subset["gap"], subset["sheafGAP"], facecolor=color, edgecolor=color, alpha=0.7, marker=".", s=1, zorder=2)
    x = subset['gap']
    y = subset['sheafGAP']
    
    if len(subset) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)
        ax.text(0.05, 0.95, f'y = {slope:.3f}x {intercept:+.3f}\n$R^2$ = {r_value**2:.3f}', transform=ax.transAxes, verticalalignment='top', fontsize=9)

    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    plt.tight_layout()
    plt.savefig(f"HASgap{atom}.pdf", format="pdf", bbox_inches='tight')
    plt.clf()
