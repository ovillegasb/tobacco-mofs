"""Tobacco Command Line Interface (CLI)."""

import argparse
from tobacco.make_topologies import download_topologies
import tobacco.tools as tools
import os


TITLE = """\033[1;36m
  _______    ____         _____
 |__   __|  |  _ \\       / ____|
    | | ___ | |_) | __ _| |     ___ ___
    | |/ _ \\|  _ < / _` | |    / __/ _ \\
    | | (_) | |_) | (_| | |___| (_| (_) |
    |_|\\___/|____/ \\__,_|\\_____\\___\\___/
\033[m
Python module tobacco or Topologically Based Crystal Constructor was developed
to rapidly produce molecular representations of porous crystals as
crystallographic information (.cif) files, which can then be used for
molecular simulation or for materials characterization.

Modifications: Orlando VILLEGAS
Date: 2025-06-03

Authors:
    Andrew S. Rosen
    Ryther Anderson
    Andrey A. Bezrukov
Date: 2021-04-08
https://github.com/tobacco-mofs/tobacco_3.0

"""

# Database path
root = os.path.dirname(__file__)
db = os.path.join(root, "data")


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
        help="Download and configure the topology database. The database\
        is obtained from 'http://rcsr.net'. You can change the default\
        url using the option `--db_url`.",
        action="store_true",
        default=False
    )

    fileinput.add_argument(
        "--db_url",
        help="Defines the url of the RCSR nets database.",
        default="http://rcsr.net/downloads/RCSRnets.cgd",
        dest="db",
        metavar="url"
    )

    fileinput.add_argument(
        "--db_path",
        help="Indicates the database path.",
        default=db,
        metavar="path"
    )

    fileinput.add_argument(
        "--check_top",
        help="Check if the topology exists.",
        default=None,
        metavar="topology"
    )

    nodeblock = parser.add_argument_group(
        "\033[1;36mGenerate node block\033[m")

    nodeblock.add_argument(
        "--show_pgs",
        help="Show available point groups.",
        action="store_true",
        default=False
    )

    nodeblock.add_argument(
        "-m", "--metal",
        help="Metalic center to cif format MXn. It must be used with the\
        -pg/--pointgroup option to choose the point group to use. Or use\
        --gen_geometries to generate a metal center with all available\
        geometries.",
        type=str,
        default=None,
        metavar="element"
    )

    nodeblock.add_argument(
        "-pg", "--pointgroup",
        help="Defines the point group used for the metallic core.",
        type=str,
        default=None,
        metavar="pointgroup"
    )

    nodeblock.add_argument(
        "--gen_geometries",
        help="Generate geometries using a metal center",
        action="store_true",
        default=False
    )

    nodeblock.add_argument(
        "-d", "--distance",
        help="Distance of the dummy atoms to the metallic center. By default\
        is 0.8 angs",
        type=float,
        default=0.8,
        metavar="r"
    )

    edgeblock = parser.add_argument_group(
        "\033[1;36mGenerate edge block\033[m")

    edgeblock.add_argument(
        "--make_XX_edge",
        help="Generate an edge without ligand, X--X. It functions as a\
        connector to more complex topologies.",
        action="store_true",
        default=False
    )

    edgeblock.add_argument(
        "-l", "--ligand",
        help="Defines the edge block from a ligand structural file. The block\
        should be accompanied with the -X/--ndx_X option indicating two atoms\
        that function as connectors.",
        type=str,
        default=None,
        metavar="structure"
    )

    edgeblock.add_argument(
        "-X", "--ndx_X",
        help="Index atoms to be used as dummy atoms.",
        type=int,
        nargs="+",
        default=None
    )

    ion_insertion = parser.add_argument_group(
        "\033[1;36mOptions used to control ion insertion.\033[m")

    ion_insertion.add_argument(
        "--ion",
        help="Ion file name.",
        type=str,
        default=None,
        metavar="file"
    )

    ion_insertion.add_argument(
        "--n_ions",
        help="Number of ions to insert.",
        type=int,
        default=0,
        metavar="N"
    )

    fileoutput = parser.add_argument_group(
        "\033[1;36mOutput settings\033[m")

    fileoutput.add_argument(
        "-o", "--output",
        help="Output file name.",
        type=str,
        default=None,
        metavar="file"
    )

    RunTobacco = parser.add_argument_group(
        "\033[1;36mOptions to run ToBaCco\033[m")

    RunTobacco.add_argument(
        "--make_MOF",
        help="Generates a MOF from the nodes and edges contained in the \
working directory from a defined topology. Use --make_MOF for default \
generation, or --make_MOF 'topol' to specify a topology.",
        nargs="?",
        const=True,
        default=False
    )

    RunTobacco.add_argument(
        "-nt", "--n_node_type",
        help="Define the number of N different node types.",
        type=int,
        default=10,
        metavar="num"
    )

    RunTobacco.add_argument(
        "-bond",
        help="Connection site bond length.",
        type=float,
        default=1.0,
        metavar="r",
        dest="connection_bond"
    )

    RunTobacco.add_argument(
        "--n_max_atoms",
        help="Maximum number of atoms allowed per structure.",
        type=int,
        default=100,
        metavar="num"
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

    RunTobacco.add_argument(
        "--z_spacing",
        help="Defined spacing for Z for 2D topologies. Default is 4.0\
        anstroms.",
        type=float,
        default=4.0,
        metavar="z_distance",
        dest="desired_z_spacing"
    )

    return vars(parser.parse_args())


def main():
    """Run main function."""
    print(TITLE)
    args = options()

    if args["get_topols_db"]:
        url = args["db"]
        download_topologies(url)

    elif args["show_pgs"]:
        tools.available_pg()

    elif args["check_top"]:
        topol = args["check_top"]

        print("Loading ToBaCco database ...")
        topols_dict = tools.load_database(db_path=args["db_path"])

        if topol in topols_dict:
            print(f"The '{topol}' topology \033[1;36mexists\033[m in the \
database.")
        else:
            print(f"The '{topol}' topology does not exist in the database.")

    elif args["metal"] is not None:
        print("Generating a node block")
        if args["pointgroup"] is not None:
            tools.gen_sbu_metal_center(**args)
        elif args["gen_geometries"] is not None:
            tools.gen_geometries_metal(**args)
        else:
            raise ValueError("Nothing has been generated because the metal\
 option comes with the point group option: use -pg/--pointgroup or use\
 --gen_geometries.")

    elif args["ligand"] is not None:
        print("Generating a edge block")
        tools.gen_sbu_edge(**args)

    elif args["make_MOF"]:
        topol = args["make_MOF"]
        db_path = args["db_path"]
        topols_dict = tools.load_database(db_path)
        templates_path = "" if db_path != db else db
        args["templates_path"] = templates_path

        print("Running ToBaCco...")

        if isinstance(topol, str):
            print("Generation using a single topology")
            print("Topology selected:", topol)
            tools.make_MOF(
                topols_dict[topol],
                n_node_type=args["n_node_type"],
                n_max_atoms=args["n_max_atoms"],
                connection_bond=args["connection_bond"],
                desired_z_spacing=args["desired_z_spacing"],
                templates_path=templates_path,
                ion=args["ion"],
                n_ions=args["n_ions"]
            )

        elif topol is True:
            print("Generation using all database topologies")
            if args["run_parallel"]:
                tools.run_tobacco_parallel(
                    topols_dict,
                    **args
                )
            else:
                tools.run_tobacco_serial(
                    topols_dict,
                    **args
                )

    return "Done!"
