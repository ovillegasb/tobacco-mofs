#!/bin/env python
# -*- coding: utf-8 -*-

"""Topologically Based Crystal Constructor (ToBaCco).

ToBaCco or Topologically Based Crystal Constructor was developed to rapidly produce molecular
representations of porous crystals as crystallographic information (.cif) files, which can then be
used for molecular simulation or for materials characterization.

Original version: https://github.com/tobacco-mofs/tobacco_3.0

This is a version of ToBaCco_3.0 intended to be used as a module, with this is looking to couple
some of its functions with other code.
"""

import os
import argparse
from tobacco.tools import build_sbu_metal_center, extract_info, write_cif_SBU
from tobacco.tools import make_MOF, load_database, build_sbu_from_gaus, run_tobacco_parallel, run_tobacco_serial
from tobacco.tools import gen_geometries_metal
from tobacco.make_topologies import topols_config, get_topols_db
from tobacco.tools import skeleton_X, make_XX_edge


TITLE = """\033[1;36m
  _______    ____         _____
 |__   __|  |  _ \\       / ____|
    | | ___ | |_) | __ _| |     ___ ___
    | |/ _ \\|  _ < / _` | |    / __/ _ \\
    | | (_) | |_) | (_| | |___| (_| (_) |
    |_|\\___/|____/ \\__,_|\\_____\\___\\___/
\033[m
Python module tobacco or Topologically Based Crystal Constructor was developed to rapidly
 produce molecular representations of porous crystals as crystallographic information (.cif) files,
 which can then be used for molecular simulation or for materials characterization.

Modifications: Orlando VILLEGAS
Date: 2024-10-15

Authors:
    Andrew S. Rosen
    Ryther Anderson
    Andrey A. Bezrukov
Date: 2021-04-08
https://github.com/tobacco-mofs/tobacco_3.0

Usage:

    Download and setup topology database:
        tobacco --get_topols_db

    To generate a cif file where a block is built with the metal in the center
    with dummy atoms forming the indicated geometry:
        tobacco -m Sc -pg Oh -d 0.8 -o 6X_Sc.cif

    Show available punctual groups.
        tobacco --show_pgs

    Generate all available geometries for a metal core:
        tobacco --gen_geometries Sc -d 1.0

    Convert .com file to SBU ToBaCco, used to create a structure with dummy atoms to be removed (X-->Fr):
        tobacco --build_sbu node -i 4X_C2.com -o 4X_C2.cif

    Method to generate cif files used for edges, used to preserve the dummy atom:
        tobacco --build_sbu edge -i 2X_SCN.com -X 1 2 -o 2X_SCN.cif

    Generate an edge without ligand, X--X. It functions as a connector to more complex topologies.
        tobacco --make_XX_edge

    Check if the topology exists:
        tobacco --check_top -t name

        Example names (n vertex):
            sql: 1
            ttth: 3

    To generate MOF using ToBaCco:
        tobacco --make_MOF -t pcu

    Run using all topologies en Parallel
        tobacco --make_MOF --all_topols --run_parallel

"""


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="tobacco",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        "\033[1;36mInitial settings\033[m")

    fileinput.add_argument(
        "--get_topols_db",
        help="Download and configure the topology database.",
        action="store_true",
        default=False
    )

    fileinput.add_argument(
        "-db", "--database",
        dest="db",
        help="RCSR topology database url address.\
(by default: http://rcsr.net/downloads/RCSRnets.cgd)",
        type=str,
        default="http://rcsr.net/downloads/RCSRnets.cgd",
        metavar="url"
    )

    fileinput.add_argument(
        "-m", "--metal",
        help="Metalic center to cif format MXn.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "--show_pgs",
        help="Show available point groups.",
        action="store_true",
        default=False
    )

    fileinput.add_argument(
        "-pg", "--pointgroup",
        help="Defines the point group used for the metallic core.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "-d", "--distance",
        help="Distance of the dummy atoms to the metallic center.",
        type=float,
        default=0.8
    )

    fileinput.add_argument(
        "--gen_geometries",
        help="Generate geometries using a metal center",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "--make_XX_edge",
        help="Generate an edge without ligand, X--X. It functions as a connector to more complex\
topologies.",
        action="store_true",
        default=False
    )

    fileinput.add_argument(
        "-i", "--input",
        help="Input file path.",
        type=str,
    )

    fileinput.add_argument(
        "-X", "--ndx_X",
        help="Index atoms to be used as dummy atoms.",
        type=int,
        nargs="+",
        default=None
    )

    fileinput.add_argument(
        "--build_sbu",
        help="Activates the secondary building unit (SBU) generation mode, the options are [node] or [edge].",
        type=str,
        default=None,
        dest="sbu",
        metavar="sbu"
    )

    fileoutput = parser.add_argument_group(
        "\033[1;36mOutput settings\033[m")

    fileoutput.add_argument(
        "-o", "--output",
        help="Output file name.",
        type=str,
        default="unknown.cif"
    )

    RunTobacco = parser.add_argument_group(
        "\033[1;36mOptions to run ToBaCco\033[m")

    RunTobacco.add_argument(
        "--make_MOF",
        help="Generates a MOF from the nodes and edges contained in the working directory from a\
defined topology.",
        action="store_true"
    )

    RunTobacco.add_argument(
        "-t", "--topology",
        help="Topology selected to generate a MOF.",
        type=str,
        default=None
    )

    RunTobacco.add_argument(
        "-nt", "--n_node_type",
        help="Define the number of N different node types.",
        type=int,
        default=10
    )

    RunTobacco.add_argument(
        "--n_max_atoms",
        help="Maximum number of atoms allowed per structure.",
        type=int,
        default=100
    )

    RunTobacco.add_argument(
        "--check_top",
        help="Check if the topology exists.",
        action="store_true"
    )

    RunTobacco.add_argument(
        "--all_topols",
        help="Run ToBaCco using all database topologies.",
        action="store_true"
    )

    RunTobacco.add_argument(
        "--run_parallel",
        help="Run ToBaCco in Parallel.",
        action="store_true"
    )

    return vars(parser.parse_args())


def main():
    """Run main function."""
    print(TITLE)
    args = options()

    if args["get_topols_db"]:
        url = args["db"]
        get_topols_db(url)
        topols_config(consider_2D=True)

    elif args["show_pgs"]:
        print("Available point groups:")
        for pg in skeleton_X:
            print(f"    {pg}")

    elif args["make_XX_edge"]:
        print("Making XX edge block:")
        if not os.path.exists("./edges"):
            os.mkdir("./edges")
        make_XX_edge()

    elif args["metal"] is not None:
        # SBU metal center
        element = args["metal"]
        pointgroup = args["pointgroup"]
        d = args["distance"]
        out_file = args["output"]
        sbu = build_sbu_metal_center(element, pointgroup, d)
        parameters = extract_info(sbu)

        if not os.path.exists("./nodes"):
            os.mkdir("./nodes")

        if not out_file.startswith("nodes"):
            out_file = os.path.join("nodes", out_file)

        # write cif
        write_cif_SBU(out_file, *parameters)

    elif args["gen_geometries"] is not None:
        if not os.path.exists("./nodes"):
            os.mkdir("./nodes")

        element = args["gen_geometries"]
        gen_geometries_metal(element, args["distance"])

    elif args["sbu"] is not None:
        assert args["sbu"] == "node" or args["sbu"] == "edge", "SBU type not selected"
        assert args["input"] is not None, "Input file not selected"
        print("A secondary building unit will be generated")
        in_file = args["input"]
        out_file = args["output"]
        ext = in_file.split(".")[-1]
        sbu_type = args["sbu"] + "s"
        ndx_X = args["ndx_X"]
        print(f"File: {in_file}, ext: {ext}, sbu type: {sbu_type}")
        if not os.path.exists(sbu_type):
            os.mkdir(sbu_type)

        # SBU type node
        if not out_file.startswith(sbu_type):
            out_file = os.path.join(sbu_type, out_file)
        print("Out file:", out_file)
        if in_file.endswith(".com"):
            sbu = build_sbu_from_gaus(in_file)
            # It checks if the file contains dummy atoms by default.
            # If there are no atoms, check if indices have been defined to indicate which atoms
            # will function as non-removable dummy atoms.
            if "X" not in sbu.get_chemical_symbols():
                if ndx_X is None:
                    print("No dummy atoms are found in the structure and no defined indices\
(non-removable X) are found.")
                    print("="*10)
                    for at in sbu:
                        print("    ", at.index, at.symbol)
                    print("="*10)
                    print("Example use: -X 0 1")
                    raise ValueError("No dummy atoms have been defined.")

                parameters = extract_info(sbu, ndx_X, remove_dummy=False)
            else:
                parameters = extract_info(sbu, remove_dummy=True)

            # write cif
            write_cif_SBU(out_file, *parameters)

    elif args["make_MOF"]:
        if not args["all_topols"]:
            topol = args["topology"]
            assert topol is not None, "Define a topology using -t option"

            print("Running ToBaCco")
            topols_dict = load_database()
            print("Topology selected:", topol)
            make_MOF(topols_dict[topol], n_node_type=args["n_node_type"], n_max_atoms=args["n_max_atoms"])

        else:
            print("Running ToBaCco")
            print("Using all topologies")
            topols_dict = load_database()
            if args["run_parallel"]:
                run_tobacco_parallel(topols_dict, n_node_type=args["n_node_type"], n_max_atoms=args["n_max_atoms"])
            else:
                run_tobacco_serial(topols_dict, n_node_type=args["n_node_type"], n_max_atoms=args["n_max_atoms"])

    elif args["check_top"]:
        topol = args["topology"]
        assert topol is not None, "Define a topology using -t option"

        print("Loading ToBaCco database ...")
        topols_dict = load_database()

        if topol in topols_dict:
            print(f"The '{topol}' topology \033[1;36mexists\033[m in the database.")
        else:
            print(f"The '{topol}' topology does not exist in the database.")


if __name__ == '__main__':
    main()
