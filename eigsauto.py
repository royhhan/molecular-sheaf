from collections import defaultdict
import numpy as np
import networkx as nx
from rdkit import Chem


smiles = input("SMILES: ")
electrons = input("Electron count: ")

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

print(eigenvalues)

degeneracies = defaultdict(int)
for val in eigenvalues:
    degeneracies[val] += 1

for val in degeneracies:
    count = degeneracies[val] * 2
    if electrons - count < 0:
        if electrons - count == -1:
            print("unstable doublet")
        elif electrons - count == -2:
            print("unstable triplet")
        else:
            print("unstable unclassed")
        break
    elif electrons - count == 0:
        print("stable singlet")
        break
    else:
        electrons -= count
