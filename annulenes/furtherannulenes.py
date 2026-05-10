import numpy as np
import networkx as nx
from rdkit import Chem
from matplotlib import pyplot as plt

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

for i in range(4, 50, 2):
    smiles = "C1" + ((i-2)//2)*"=CC" + "=C1"
    electrons = i

    eigenvalues = findEigs(smiles)

    homo = 0
    lumo = 0
    for val in eigenvalues:

        if homo > 0:
            lumo = val
            x.append(i)
            y.append(lumo-homo)
            break
        elif electrons - 2 <= 0:
            homo = val
        else:
            electrons -= 2
            
plt.figure(figsize=(10,6))

plt.xlabel('Number of carbons in annulene')
plt.ylabel('Frontier eigengap')
plt.grid(visible=True, which="both", color="0.8", zorder=1)
plt.scatter(x, y, facecolor=(0.0078, 0.8000, 0.9961, 0.5), edgecolor="xkcd:bright sky blue", zorder=2)
for i in range(len(x)):
    plt.annotate(x[i], xy=(x[i], y[i]), xytext=(5, 5), textcoords='offset points', color="0", size="medium", ha="left")

plt.xlim(left=0)
plt.ylim(bottom=0)
plt.savefig("furtherannulenes.pdf", format="pdf", bbox_inches='tight')
plt.clf()
