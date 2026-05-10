import numpy as np
import networkx as nx
from rdkit import Chem
import pandas as pd

def findEigs(mol):
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


def sheafCompute(smiles):
    mol = Chem.MolFromSmiles(smiles)

    electrons = 0
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        
        if sym == 'C': 
            electrons += 1
        elif sym == 'B': 
            electrons += 0
        elif sym in ['O', 'S']: 
            electrons += 2
        elif sym == 'N':
            if atom.GetDegree() == 3 or atom.GetTotalNumHs() > 0:
                electrons += 2
            else:
                electrons += 1
        else:
            print("unrecognised")
            electrons += 1 

    n = electrons

    eigenvalues = findEigs(mol)

    homo = 0
    lumo = 0
    for val in eigenvalues:

        if homo > 0:
            lumo = val
            break
        elif electrons - 2 <= 0:
            homo = val
        else:
            electrons -= 2

    return pd.Series({
        "sheafHOMO": homo,
        "sheafLUMO": lumo,
        "sheafGAP": lumo-homo,
        "electrons": n
    })


# COMPAS datasets sourced from https://gitlab.com/porannegroup/compas
allpass = pd.read_csv("compas-2D.csv", usecols=["smiles", "charge", "homo", "lumo", "gap", "b", "s", "o", "n"])
allpass = allpass[allpass["charge"] == 0]

allpass[["sheafHOMO", "sheafLUMO", "sheafGAP", "electrons"]] = allpass['smiles'].apply(sheafCompute)

allpass.to_csv("sheafhass.csv", index=False)

mask = (allpass[["b", "s", "o", "n"]] != 0).sum(axis=1) == 1
onehetero = allpass[mask]
onehetero.to_csv("sheafhasshetero.csv", index=False)
