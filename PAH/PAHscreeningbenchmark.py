import time
import tracemalloc
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem

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


def benchmark_molecule(smiles, repeats=10):
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        sheafCompute(smiles)
        times.append(time.perf_counter() - t0)

    tracemalloc.start()
    sheafCompute(smiles)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "mean_s":   np.mean(times),
        "std_s":    np.std(times),
        "peak_kb":  peak / 1024
    }

def run_full_benchmark(df, repeats=5):
    records = []
    for _, row in df.iterrows():
        smiles = row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        n_atoms = mol.GetNumAtoms()
        n_bonds = mol.GetNumBonds()

        stats = benchmark_molecule(smiles, repeats=repeats)
        records.append({
            "smiles":   smiles,
            "n_atoms":  n_atoms,
            "n_bonds":  n_bonds,
            **stats
        })

    return pd.DataFrame(records)



cata = pd.read_csv("compas-1D_cam-b3lyp_aug-cc-pvdz.csv", usecols=["smiles", "HOMO_eV", "LUMO_eV", "GAP_eV"])
peri = pd.read_csv("compas-3D.csv", usecols=["smiles", "HOMO_eV", "LUMO_eV", "GAP_eV"])
allpahs = pd.concat([cata, peri], ignore_index=True)
bdf = run_full_benchmark(allpahs, repeats=5)

print(bdf[["n_atoms", "mean_s", "std_s", "peak_kb"]].describe().round(4))
