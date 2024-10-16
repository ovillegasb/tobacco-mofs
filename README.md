
# TOBACCO modified version

*Orlando Villegas* - **2024**

Python module tobacco or Topologically Based Crystal Constructor was developed to rapidly produce molecular representations of porous crystals as crystallographic information (.cif) files, which can then be used for molecular simulation or for materials characterization.

This version has as starting point the tobacco_3.0 (https://github.com/tobacco-mofs/tobacco_3.0) version, I have made different modifications to make it executable from command line.

This version keeps the same license (General Public License v3.0) and is modifiable. What I have done is to adapt it to my PostDoc work and share it with the community.

The main use I gave it was as a topological MOF generator. Topological MOF generator by molecular block alignment (ToMOFGen)

## Installation

You can start from a virtual environment local:

```
python -m venv .venv --copies --prompt Tobacco
source .venv/bin/activate
pip install ToBaCco-Modified
```

## Pre-use

Show installation information:

```
pip show tobacco
```

To check that tobacco is installed correctly you can run:

```
tobacco -h
```

Before you can use tobacco you must download the topology database. To do this, execute, this may take a few minutes:

```
tobacco --get_topols_db
```

## Example of uses

To generate a cif file where a block is built with the metal in the center with dummy atoms forming the indicated geometry:

    tobacco -m Sc -pg Oh -d 1.0 -o 6X_Sc.cif


Show available punctual groups.

    tobacco --show_pgs


Generate all available geometries for a metal core:

    tobacco --gen_geometries Sc -d 1.0


Convert `.com` file (Gaussian format) to SBU (Secundary Building Unit) ToBaCco, used to create a structure with dummy atoms to be removed (X-->Fr):

    tobacco --build_sbu node -i 4X_C2.com -o 4X_C2.cif


Method to generate cif files used for edges, and used to preserve the dummy atom:

    tobacco --build_sbu edge -i 2X_SCN.com -X 1 2 -o 2X_SCN.cif


Generate an edge without ligand, X--X. It functions as a connector to more complex topologies.

    tobacco --make_XX_edge


Check if the topology exists:

    tobacco --check_top -t pcu


To generate MOF using ToBaCco:

    tobacco --make_MOF -t pcu


Run using all topologies en Parallel

    tobacco --make_MOF --all_topols --run_parallel

## Remove

For remove tobacco, use:

```
pip uninstall tobacco
```