# ToBaCCo-MOFS

![version](https://img.shields.io/badge/version-3.1.2-blue)


Topologically Based Crystal Constructor (**ToBaCCo**) is a Python tool developed to rapidly generate molecular representations of porous crystals as Crystallographic Information Files (`.cif`). These can be used for molecular simulation and materials characterization.

## Version

This is version **3.1.2**, a modified and extended version of [ToBaCCo_3.0](https://github.com/tobacco-mofs/tobacco_3.0).

- Original Authors:
  - Ryther Anderson
  - Yamil Colón
  - Diego Gómez-Gualdrón

- Modified and extended by:
  - **Orlando Villegas** (2025)


## Key Features in v3.1.2

- Command-line interface (CLI) for simplified use.
- Extended topology generation options.
- New routines for SBU generation from Gaussian `.com` files.
- Integrated topology database handling.
- Better Python 3 compatibility.

---


## Installation

You can install it from source or package (if uploaded to PyPI):

### Virtual Environment Setup (recommended)

```bash
python -m venv .venv --copies --prompt tobacco
source .venv/bin/activate
```

## Getting Started


### Check installation:


```bash
tobacco -h
```

### Download the topology database:

Before you can use tobacco you must download the topology database. To do this, execute, this may take a few minutes:

```bash
tobacco --get_topols_db
```

## Example of uses

1. Generate geometry around a metal center:

To generate a cif file where a block is built with the metal in the center with dummy atoms forming the indicated geometry:

```bash
tobacco -m Sc -pg Oh -d 1.0 -o 6X_Sc.cif
```

2. Show available point groups:

```bash
tobacco --show_pgs
```

3. Generate geometries for a metal:

```bash
tobacco --gen_geometries Sc -d 1.0
```

4. Build SBU from `.com` file (Gaussian):

Convert `.com` file (Gaussian format) to SBU (Secundary Building Unit) ToBaCco, used to create a structure with dummy atoms to be removed (X-->Fr):

```bash
tobacco --build_sbu node -i 4X_C2.com -o 4X_C2.cif
```

5. Build edge SBU:

Method to generate cif files used for edges, and used to preserve the dummy atom:

```bash
tobacco --build_sbu edge -i 2X_SCN.com -X 1 2 -o 2X_SCN.cif
```


6. Generate a dummy edge (X--X):

Generate an edge without ligand, X--X. It functions as a connector to more complex topologies.

```bash
tobacco --make_XX_edge
```


7. Check topology availability:

```bash
tobacco --check_top -t pcu
```

8. Generate MOF structure with given topology:


```bash
tobacco --make_MOF -t pcu
```

9. Generate MOFs for all topologies in parallel:


```bash
tobacco --make_MOF --all_topols --run_parallel
```


## 📄 License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](./LICENSE) file for more details.

This version is a fork and modification of ToBaCCo_3.0, developed by Ryther Anderson et al.