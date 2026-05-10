import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats


allpahs = pd.read_csv("sheafpahs.csv")


# gap
plt.figure(figsize=(10,6))

plt.xlabel('DFT HOMO-LUMO gap / eV')
plt.ylabel('Frontier eigengap')
plt.grid(visible=True, which="both", color="0.8", zorder=1)
plt.scatter(allpahs["GAP_eV"], allpahs["sheafGAP"], facecolor=(0.0078, 0.8000, 0.9961, 0.5), edgecolor="xkcd:bright sky blue", alpha=0.3, marker=".", s=1, zorder=2)

x = allpahs['GAP_eV']
y = allpahs['sheafGAP']

# compute line of best fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

# plot line of best fit
x_line = np.linspace(x.min(), x.max(), 100)
plt.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)

# add equation and r^2
plt.text(0.05, 0.95, f'y = {slope:.3f}x + {intercept:.3f}\n$R^2$ = {r_value**2:.3f}', transform=plt.gca().transAxes, verticalalignment='top', fontsize=10)

plt.savefig("PAHgap.pdf", format="pdf", bbox_inches='tight')
plt.clf()

# density map 
plt.figure(figsize=(10,6))

plt.xlabel('DFT HOMO-LUMO gap / eV')
plt.ylabel('Frontier eigengap')
plt.grid(visible=True, which="both", color="0.8", zorder=1)

x = allpahs['GAP_eV']
y = allpahs['sheafGAP']

cmap = plt.cm.viridis_r.copy()
cmap.set_under('none')
plt.hist2d(allpahs['GAP_eV'], allpahs['sheafGAP'], bins=200, cmap=cmap, vmin=1, zorder=2)
plt.colorbar(label='Count')

# compute line of best fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

# plot line of best fit
x_line = np.linspace(x.min(), x.max(), 100)
plt.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)

# add equation and r^2
plt.text(0.05, 0.95, f'y = {slope:.3f}x + {intercept:.3f}\n$R^2$ = {r_value**2:.3f}', transform=plt.gca().transAxes, verticalalignment='top', fontsize=10)

plt.savefig("PAHgapdensitysquares.pdf", format="pdf", bbox_inches='tight')
plt.clf()
