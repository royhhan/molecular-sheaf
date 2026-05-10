Code and data for **A Cellular Sheaf Model for Aromaticity in Molecular Systems**, submitted to the Rouse Research Essay Award 2026. This repository contains all scripts needed to reproduce the results in the paper.

---

## Repository structure

```
molecular-sheaf/
├── README.md
├── coboundaries.py              # Coboundaries of benzene & cyclobutadiene
├── eigsauto.py                  # SMILES → eigenvalues computation
├── eigsmanual.py                # Manual eigenvalue computation
├── annulenes/
│   ├── furtherannulenes.py      # Frontier eigengap across annulenes [4]–[48]
│   ├── furtherannulenes.pdf
│   ├── furtherannulenes2.py     # Frontier eigengap vs ISE for select annulenes
│   └── furtherannulenes2.pdf
├── PAH/
│   ├── compas-1D_cam-b3lyp_aug-cc-pvdz.csv
│   ├── compas-3D.csv
│   ├── sheafpahs.csv            # Precomputed sheaf results
│   ├── PAHscreening.py          # Sheaf computation from COMPAS-1D and COMPAS-3D
│   ├── PAHscreeningbenchmark.py 
│   ├── PAHgraphing.py           # Scatter and density plots of ΔF vs DFT gap
│   ├── PAHgap.pdf
│   ├── PAHgapdensitysquares.pdf
│   └── PAHgraphinginteractive.py
└── HAS/
    ├── compas-2D.csv
    ├── sheafhass.csv            # Precomputed sheaf results
    ├── sheafhasshetero.csv      # Precomputed sheaf results, single heteroatoms
    ├── HASscreening.py          # Sheaf computation from COMPAS-2D
    ├── HASscreeningbenchmark.py
    ├── HASgraphing.py           # Scatter and density plots of ΔF vs DFT gap
    ├── HASgap.pdf
    ├── HASgapdensitysquares.pdf
    ├── HASgraphingbyheteroatom.py
    ├── HASgaphetero.pdf         # All heteroatoms combined, coloured by type
    ├── HASgapb.pdf              # Boron-containing molecules
    ├── HASgapn.pdf              # Nitrogen-containing molecules
    ├── HASgapo.pdf              # Oxygen-containing molecules
    └── HASgaps.pdf              # Sulfur-containing molecules
```

---

## Reproducing the results

### Dependencies

```
numpy
scipy
matplotlib
networkx
rdkit
pandas
```

### Model systems (Section 3.1)

`coboundaries.py` gives the manually constructed coboundary matrices for benzene and cyclobutadiene, which can then be processed in `eigsmanual.py`. `eigsauto.py`accepts any SMILES string interactively and returns the sheaf eigenvalues and electronic stability classification.

### Annulene series (Section 3.1)

Run from the `annulenes/` directory:

```
python furtherannulenes.py
python furtherannulenes2.py
```

### PAH large-scale validation (Section 3.2)

Run from the `PAH/` directory. `sheafpahs.csv` is precomputed using `PAHscreening.py` and can be used directly with the graphing scripts. To recompute from scratch:

```
python PAHscreening.py
python PAHgraphing.py
```

### Heteroaromatic screening (Section 3.3)

Run from the `HAS/` directory. `sheafhass.csv` and `sheafhasshetero.csv` are precomputed using `HASscreening.py` and can be used directly with the graphing scripts. To recompute from scratch:

```
python HASscreening.py
python HASgraphing.py
python HASgraphingbyheteroatom.py
```

---

## Data

Annulene ISEs are sourced from https://pubs.acs.org/doi/10.1021/ol027571b.

Raw COMPAS data is sourced from the [COMPAS Project](https://gitlab.com/porannegroup/compas):

- **COMPAS-1D** — cata-condensed polybenzenoid hydrocarbons
- **COMPAS-2D** — heteroatom-substituted polycyclic aromatic systems
- **COMPAS-3D** — peri-condensed polybenzenoid hydrocarbons
