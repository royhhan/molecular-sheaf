import numpy as np
import networkx as nx
from rdkit import Chem
import pandas as pd

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


def sheafCompute(smiles):
    mol = Chem.MolFromSmiles(smiles)
    electrons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    n = electrons

    eigenvalues = findEigs(smiles)

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
cata = pd.read_csv("compas-1D_cam-b3lyp_aug-cc-pvdz.csv", usecols=["smiles", "HOMO_eV", "LUMO_eV", "GAP_eV"])
peri = pd.read_csv("compas-3D.csv", usecols=["smiles", "HOMO_eV", "LUMO_eV", "GAP_eV"])
allpahs = pd.concat([cata, peri], ignore_index=True)

allpahs[["sheafHOMO", "sheafLUMO", "sheafGAP", "electrons"]] = allpahs['smiles'].apply(sheafCompute)

allpahs.to_csv("sheafpahs.csv", index=False)
