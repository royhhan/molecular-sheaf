import time
import tracemalloc
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem

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


allpass = pd.read_csv("compas-2D.csv", usecols=["smiles", "charge", "homo", "lumo", "gap", "b", "s", "o", "n"])
allpass = allpass[allpass["charge"] == 0]
bdf = run_full_benchmark(allpass, repeats=5)

print(bdf[["n_atoms", "mean_s", "std_s", "peak_kb"]].describe().round(4))
