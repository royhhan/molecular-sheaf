import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats

allpahs = pd.read_csv("sheafpahs.csv")

# gap
fig, ax = plt.subplots(figsize=(10,6))
ax.set_xlabel('DFT HOMO-LUMO gap / eV')
ax.set_ylabel('Frontier eigengap')
ax.grid(visible=True, which="both", color="0.8", zorder=1)
sc = ax.scatter(allpahs["GAP_eV"], allpahs["sheafGAP"], facecolor=(0.0078, 0.8000, 0.9961, 0.5), edgecolor="xkcd:bright sky blue", alpha=0.3, marker=".", s=1, zorder=2)
x = allpahs['GAP_eV']
y = allpahs['sheafGAP']
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
x_line = np.linspace(x.min(), x.max(), 100)
ax.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)
ax.text(0.05, 0.95,
        f'y = {slope:.3f}x + {intercept:.3f}\n$R^2$ = {r_value**2:.3f}',
        transform=ax.transAxes,
        verticalalignment='top',
        fontsize=10)

def on_click(event):
    if event.inaxes:
        distances = np.hypot(allpahs['GAP_eV'] - event.xdata,
                             allpahs['sheafGAP'] - event.ydata)
        idx = distances.idxmin()
        row = allpahs.iloc[idx]
        print(f"SMILES: {row['smiles']}")
        print(f"GAP: {row['GAP_eV']:.3f} eV")
        print(f"sheafGAP: {row['sheafGAP']:.3f}")
        print("-" * 40)

fig.canvas.mpl_connect('button_press_event', on_click)

plt.show()
plt.clf()
