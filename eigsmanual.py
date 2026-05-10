import numpy as np
import scipy.sparse as sp
from collections import defaultdict

import coboundaries

# inputs 
electrons = 30
coboundary = coboundaries.C6H6


coboundary = sp.csr_matrix(coboundary)

laplacian = coboundary.T @ coboundary
size = laplacian.shape[0]

laplacian = laplacian.toarray()
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
