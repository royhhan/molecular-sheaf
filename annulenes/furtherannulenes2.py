import numpy as np
import networkx as nx
from rdkit import Chem
from matplotlib import pyplot as plt
from scipy import stats


def findEigs(smiles):
    mol = Chem.MolFromSmiles(smiles)

    adjacency = np.array(Chem.rdmolops.GetAdjacencyMatrix(mol))

    graph = nx.from_numpy_array(adjacency)
    inc_matrix_sparse = nx.incidence_matrix(graph, oriented=True)
    coboundary = inc_matrix_sparse.T.todense()

    laplacian = coboundary.T @ coboundary

    eigenvalues = np.linalg.eigvalsh(laplacian)

    eigenvalues = eigenvalues.tolist()
    eigenvalues = sorted([round(val) if close else round(val, 4) 
                        for close, val in 
                        zip(np.isclose(eigenvalues, np.round(eigenvalues)), eigenvalues)])

    return eigenvalues

x, y = [], []
allowed = [6, 10, 14, 18, 22, 26, 30, 42, 54, 66]

# ISEs sourced from https://pubs.acs.org/doi/10.1021/ol027571b
ise_per_pi = {6:  34.6/6,
             10: 32.6/10,
             14: 26.7/14,
             18: 27.4/18,
             22: 26.3/22,
             26: 24.9/26,
             30: 23.6/30,
             42: 22.9/42,
             54: 22.5/54,
             66: 22.4/66}

for i in allowed:
    smiles = "C1" + ((i-2)//2)*"=CC" + "=C1"
    electrons = i

    eigenvalues = findEigs(smiles)

    homo = 0
    lumo = 0
    for val in eigenvalues:

        if homo > 0:
            lumo = val
            x.append(ise_per_pi[i])
            y.append(lumo-homo)
            break
        elif electrons - 2 <= 0:
            homo = val
        else:
            electrons -= 2
            

plt.figure(figsize=(10,6))

plt.xlabel('ISE per π-electron (kcal/mol-electron)')
plt.ylabel('Frontier eigengap')
plt.grid(visible=True, which="both", color="0.8", zorder=1)
plt.scatter(x, y, facecolor=(0.0078, 0.8000, 0.9961, 0.5), edgecolor="xkcd:bright sky blue", zorder=2)
for i in range(len(x)):
    plt.annotate(allowed[i], xy=(x[i], y[i]), xytext=(5, -5), textcoords='offset points', color="0", size="medium", ha="left")


# compute line of best fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

# plot line of best fit
x_line = np.linspace(min(x), max(x), 100)
plt.plot(x_line, slope * x_line + intercept, color='xkcd:dark red', linewidth=1)

# add equation and r^2
plt.text(0.05, 0.95, f'y = {slope:.3f}x + {intercept:.3f}\n$R^2$ = {r_value**2:.3f}', transform=plt.gca().transAxes, verticalalignment='top', fontsize=10)

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig("furtherannulenes2.pdf", format="pdf", bbox_inches='tight')
plt.clf()
