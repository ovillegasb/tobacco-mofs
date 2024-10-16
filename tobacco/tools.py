"""
Submodulo de tobacco creado para compilar algunas funciones.

Las funciones presentes estan dirigidas para usar a tobacco como un modulo, y facilitan su
ejecucion a nivel de terminal.

Author: Orlando Villegas

"""

import os
import re
import datetime
import itertools
import glob
import shutil
import hashlib
import multiprocessing
from multiprocessing import Pool
import numpy as np
import pandas as pd
import ase
from ase.geometry.analysis import Analysis
import ase.io

# Specific function from ToBaCco
from tobacco.ciftemplate2graph import ct2g
from tobacco.bbcif_properties import cncalc, bbelems
from tobacco.vertex_edge_assign import vertex_assign, assign_node_vecs2edges
from tobacco.cycle_cocyle import cycle_cocyle, Bstar_alpha
from tobacco.SBU_geometry import SBU_coords
from tobacco.scale import scale
from tobacco.scaled_embedding2coords import omega2coords
from tobacco.place_bbs import scaled_node_and_edge_vectors, place_nodes, place_edges
from tobacco.adjust_edges import adjust_edges
from tobacco.remove_dummy_atoms import remove_Fr
from tobacco.write_cifs import bond_connected_components, fix_bond_sym, write_cif, distance_search_bond
from tobacco.remove_net_charge import fix_charges
from tobacco.tobacco import configuration
from tobacco.tobacco import metal_elements
from tobacco.tobacco import vname_dict


# Database path
root = os.path.dirname(__file__)
db = os.path.join(root, "data")

# ToBaCco configurations
USER_SPECIFIED_NODE_ASSIGNMENT = configuration.USER_SPECIFIED_NODE_ASSIGNMENT
# SYMMETRY_TOL = configuration.SYMMETRY_TOL
SYMMETRY_TOL = {
    2: 0.10,
    3: 0.50,  # 0.12
    4: 0.35,
    5: 0.25,
    6: 0.45,
    7: 0.35,
    8: 0.40,
    9: 0.60,
    10: 0.60,
    12: 0.60
}
ALL_NODE_COMBINATIONS = configuration.ALL_NODE_COMBINATIONS
COMBINATORIAL_EDGE_ASSIGNMENT = configuration.COMBINATORIAL_EDGE_ASSIGNMENT
SINGLE_METAL_MOFS_ONLY = configuration.SINGLE_METAL_MOFS_ONLY
MOFS_ONLY = configuration.MOFS_ONLY
CONNECTION_SITE_BOND_LENGTH = 1.5  # configuration.CONNECTION_SITE_BOND_LENGTH
FIX_UC = configuration.FIX_UC
SCALING_ITERATIONS = configuration.SCALING_ITERATIONS
PRE_SCALE = configuration.PRE_SCALE
MIN_CELL_LENGTH = configuration.MIN_CELL_LENGTH
OPT_METHOD = configuration.OPT_METHOD
WRITE_CHECK_FILES = configuration.WRITE_CHECK_FILES
CHARGES = configuration.CHARGES
ORIENTATION_DEPENDENT_NODES = configuration.ORIENTATION_DEPENDENT_NODES
RECORD_CALLBACK = configuration.RECORD_CALLBACK
PLACE_EDGES_BETWEEN_CONNECTION_POINTS = configuration.PLACE_EDGES_BETWEEN_CONNECTION_POINTS
REMOVE_DUMMY_ATOMS = configuration.REMOVE_DUMMY_ATOMS
BOND_TOL = 3.0  # configuration.BOND_TOL
WRITE_CIF = configuration.WRITE_CIF
IGNORE_ALL_ERRORS = configuration.IGNORE_ALL_ERRORS

pi = np.pi


# Skeleton of dummy atoms formed a geometry
skeleton_X = {
    "D*h": [
        "XX",
        [(-1, 0, 0), (1, 0, 0)]
    ],
    "D4h": [
        "XXXX",
        [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]
    ],
    "D2h": [
        "XXXX",
        [
            (0.86602540, -0.5, 0),
            (0.86602540, 0.5, 0),
            (-0.86602540, -0.5, 0),
            (-0.86602540, 0.5, 0)
        ]
    ],
    "D3h": [
        "XXX",
        [
            (0, 1, 0),
            (0, -np.sin(np.deg2rad(30.0)), -np.cos(np.deg2rad(30.0))),
            (0, -np.sin(np.deg2rad(30.0)), np.cos(np.deg2rad(30.0)))
        ]
    ],
    "Td": [
        "XXXX",
        [
            (0, 0, 1),
            (0, -0.94280904, -0.33333333),
            (0.81649658, 0.47140452, -0.33333333),
            (-0.81649658, 0.47140452, -0.33333333)
        ]
    ],
    "Oh": [
        "XXXXXX",
        [
            (0, 0, 1),
            (0, 1, 0),
            (1, 0, 0),
            (0, 0, -1),
            (0, -1, 0),
            (-1, 0, 0)
        ]
    ],
    "C2v": [
        "XX",
        [
            (0, 1, 0),
            (0.94551858, -0.32556815, 0)
        ]
    ],
    "C3v": [
        "XXX",
        [
            (0.8164967,  0.4714046,   0.3333331),
            (0.0000000,  -0.9428091,  0.3333331),
            (-0.8164967, 0.4714046,   0.3333331)
        ]
    ],
    "C2v-T-shape": [
        "XXX",
        [
            (1.0,  0.0,   0.0),
            (0.0,  0.0,  1.0),
            (-1.0,  0.0,   0.0)
        ]
    ],
    "C2": [
        "XXXX",
        [
            (0.4392342,  -0.0008507,  0.8983722),
            (-0.4392342, 0.0008507,   0.8983722),
            (-0.0000001, -0.4392350,  -0.8983722),
            (0.0000001,  0.4392350,   -0.8983722)
        ]
    ],
    "D3h-prism": [
        "XXXXXX",
        [
            (0.3396886,  -0.1961193,  0.9198635),
            (-0.3396886, -0.1961193,  0.9198635),
            (0.0000000,  0.3922386,   0.9198635),
            (0.0000000,  0.3922386,   -0.9198635),
            (0.3396886,  -0.1961193,  -0.9198635),
            (-0.3396886, -0.1961193,  -0.9198635)
        ]
    ],
    "S6": [
        "XXXXXX",
        [
            (0.5475622,  -0.6480586,  0.5293351),
            (-0.2874541, -0.7982321,  -0.5293351),
            (-0.8350163, -0.1501734,  0.5293351),
            (-0.5475622, 0.6480586,   -0.5293351),
            (0.2874541,  0.7982321,   0.5293351),
            (0.8350163,  0.1501734,   -0.5293351)
        ]
    ]
}


def load_database():
    """Load ToBaCco database."""
    print("ToBaCco path:", db)
    topols_pth = os.path.join(db, 'template_database/*.cif')
    topols_pth_2D = os.path.join(db, 'template_2D_database/*.cif')
    templates = sorted(glob.glob(topols_pth_2D) + glob.glob(topols_pth))
    # topols_dict = {cif.replace(".cif", ""): cif for cif in templates}
    topols_dict = {os.path.basename(cif).replace(".cif", ""): cif for cif in templates}
    return topols_dict


def dummy_geometry(poitgroup):
    """
    Return a dictionary with the chemical formula and its coordinates.

    This is done without occupying the center of the geometry. Each dummy atom has a distance of
    1.0 from the center.

    Parameters:
    -----------
    poitgroup : str
        Specifies the point group to use.
    """
    info = skeleton_X[poitgroup].copy()

    formula = info[0]
    coord = info[1]

    return formula, coord


def translate_center_to(coord, newcenter, box, box_init=[0, 0, 0]):
    """Translate center of a coordinate to a defined position."""
    newcoord = np.zeros(coord.shape)
    center_geom = np.sum(coord, axis=0) / len(coord)
    for i, atom in enumerate(coord):
        for j, q in enumerate(atom):
            L = abs(box[j] - box_init[j])
            box_center = (box[j] + box_init[j]) / 2
            dq = q - center_geom[j] + newcenter[j]
            if dq > L * 0.5 + box_center:
                dq -= L
            if dq <= -L * 0.5 + box_center:
                dq += L

            newcoord[i, j] = dq

    return newcoord


edge_XX_file = """data_3D Atomistic (2)
_audit_creation_date              2018-06-20
_audit_creation_method            'Materials Studio'
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                    10.0000
_cell_length_b                    10.0000
_cell_length_c                    10.0000
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
X1    Fr   -0.50817   0.13462   0.00000   0.00000  Uiso   1.00
X2    Fr   -0.50683   0.13509   0.00000   0.00000  Uiso   1.00
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
X1    X2     0.014   .     S
"""


def make_XX_edge():
    """Make X--X edge."""
    with open("./edges/ntn_edge.cif", "w") as f:
        f.write(edge_XX_file)
    print("XX edge saved!")


def build_sbu_metal_center(element, pointgroup, d=1.0, box=[10, 10, 10], angles=[90., 90., 90.], pbc=True):
    """
    Create an SBU from a metallic element and its point group.

    This function returns an ase.atoms object with the dimensions of the box already defined. The
    molecule created will be centered in the center of the box.

    Parameters:
    -----------
    element : str
        Metal element.
    pointgroup : str
        Desired point group.
    d : float
        Distance of the dummy atoms to the metalic center.
    box : list(3x1) of floats
        Box dimensions.
    angles L list(3x1) of floats
        Cell angles.
    pbc : bool
        Activate periodic conditions.
    """
    formula, pos = dummy_geometry(pointgroup)
    sbu = ase.Atoms(formula, positions=pos)
    try:
        sbu += ase.Atom(element, (0, 0, 0))
    except KeyError:
        pass

    sbu.positions *= d

    # box parameters
    sbu.set_cell(box + angles)
    sbu.set_pbc(pbc)

    box = np.array(box)
    sbu.positions = translate_center_to(sbu.positions, box/2, box)

    return sbu


def build_sbu_edge(file, box=[10, 10, 10], angles=[90., 90., 90.], pbc=True):
    """TODO"""
    mol = ase.io.read(file)
    mol.set_cell(box + angles)
    mol.set_pbc(pbc)

    box = np.array(box)
    mol.positions = translate_center_to(mol.positions, box/2, box)

    return mol


def build_sbu_from_gaus(file, box=[10, 10, 10], angles=[90., 90., 90.], pbc=True):
    """Construct ASE Atoms object from a Gaussian .com input file."""
    dfatoms = read_input_gaus(file)
    mol = ase.Atoms()
    for i in dfatoms.index:
        at = dfatoms.loc[i, "atsb"]
        pos = dfatoms.loc[i, ["x", "y", "z"]].values
        mol += ase.Atom(at, position=pos)

    mol.set_cell(box + angles)
    mol.set_pbc(pbc)
    box = np.array(box)
    mol.positions = translate_center_to(mol.positions, box/2, box)

    return mol


def extract_info(atoms, ndx_X=[], remove_dummy=True):
    """
    Extract information from a structure ASE.Atoms.

    Parameters:
    -----------
    atoms : ASE.Atoms
        Atoms structure of the ASE module.

    ndx_X : list
        List defining the index of atoms to be used as dummy atoms X, however, these will be
        preserved in the final mof.

    remove_dummy : bool
        Substitute the atomic symbol for Fr, which is used by ToBaCco to remove these dummy atoms.

    """
    print("Chemical formula:", atoms.symbols.get_chemical_formula(), "\n")

    # Create dict of index: X
    ndx_X = {i: "X" for i in ndx_X}

    # scaled_params
    box_params = list(atoms.cell.cellpar())
    # print("box_params:", box_params, "\n")

    # sc_unit_cell
    cell = atoms.cell.array
    # print("cell:\n", cell, "\n")

    # placed_all
    coord_list = []
    for i, at in enumerate(atoms):
        try:
            symb = ndx_X[i]
        except KeyError:
            symb = at.symbol

        if symb == "X" and remove_dummy:
            coord_list.append([f"{symb}{i+1}", "Fr", at.x, at.y, at.z])
        else:
            coord_list.append([f"{symb}{i+1}", at.symbol, at.x, at.y, at.z])

    # print("coord_list:")
    # for row in coord_list:
    #     print(row)
    # print()

    # fixed_bonds
    fixed_bonds = []
    matrixD = atoms.get_all_distances()
    ana = Analysis(atoms)
    bonds = ana.all_bonds
    for i, b in enumerate(bonds[0]):
        for ib in b:
            if "X" in coord_list[i][0] and "X" in coord_list[ib][0]:
                continue
            else:
                fixed_bonds.append([coord_list[i][0], coord_list[ib][0], matrixD[i][ib], ".", "S"])

    # print("fixed_bonds:")
    # for row in fixed_bonds:
    #     print(row)
    # print()

    return coord_list, fixed_bonds, box_params, cell


def write_cif_SBU(file, coord_list, fixed_bonds, box_params, cell, wrap_coords=True):
    """
    Save cif file.

    This cif save file uses the cif file format implemented by tobacco.

    Parameters:
    -----------
    file : str
        Output file name.

    coord_list : list(Nx5)
        Type of atom, atomic symbol and coordinates in space.

    fixed_bonds : list(Nx5)
        Type of atoms involved in the bond, bond distance, geometry and bond order.

    box_params : list(1x6) of float
        Contains the parameters of the box: a, b, c, alpha, beta, gamma.

    cell : numpy.array(3x3)
        Cell Matrix.

    """
    sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma = box_params

    lines = ""
    lines += 'data_' + "cifname" + '\n'  # + cifname[0:-4]
    lines += '_audit_creation_date              ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n'
    lines += "_audit_creation_method            'tobacco_3.0'" + '\n'
    lines += "_symmetry_space_group_name_H-M    'P1'" + '\n'
    lines += '_symmetry_Int_Tables_number       1' + '\n'
    lines += '_symmetry_cell_setting            triclinic' + '\n'
    lines += 'loop_' + '\n'
    lines += '_symmetry_equiv_pos_as_xyz' + '\n'
    lines += '  x,y,z' + '\n'
    lines += '_cell_length_a                    ' + str(sc_a) + '\n'
    lines += '_cell_length_b                    ' + str(sc_b) + '\n'
    lines += '_cell_length_c                    ' + str(sc_c) + '\n'
    lines += '_cell_angle_alpha                 ' + str(sc_alpha) + '\n'
    lines += '_cell_angle_beta                  ' + str(sc_beta) + '\n'
    lines += '_cell_angle_gamma                 ' + str(sc_gamma) + '\n'
    lines += 'loop_' + '\n'
    lines += '_atom_site_label' + '\n'
    lines += '_atom_site_type_symbol' + '\n'
    lines += '_atom_site_fract_x' + '\n'
    lines += '_atom_site_fract_y' + '\n'
    lines += '_atom_site_fract_z' + '\n'

    for row in coord_list:
        vec = list(map(float, row[2:5]))
        cvec = np.dot(np.linalg.inv(cell), vec)

        if wrap_coords:
            cvec = np.mod(cvec, 1)  # makes sure that all fractional coordinates are in [0,1]

        lines += '{:7} {:>4} {:>15} {:>15} {:>15}\n'.format(
            row[0],
            re.sub('[0-9]', '', row[1]),
            "%.10f" % np.round(cvec[0], 10),
            "%.10f" % np.round(cvec[1], 10),
            "%.10f" % np.round(cvec[2], 10)
        )

    lines += 'loop_' + '\n'
    lines += '_geom_bond_atom_site_label_1' + '\n'
    lines += '_geom_bond_atom_site_label_2' + '\n'
    lines += '_geom_bond_distance' + '\n'
    lines += '_geom_bond_site_symmetry_2' + '\n'
    lines += '_ccdc_geom_bond_type' + '\n'

    for bond in fixed_bonds:
        lines += '{:7} {:>7} {:>5} {:>7} {:>3}\n'.format(
            bond[0],
            bond[1],
            "%.3f" % float(bond[2]),
            bond[3],
            bond[4]
        )

    with open(file, "w") as CIF:
        CIF.write(lines)

    print(f"File written: {file}")


def save_state(coord_list, box, angles, name="test", ft="cif", pbc=True):
    """Save structure o cif file sing ASE."""
    struct = ase.Atoms()
    for pn in coord_list:
        at = re.sub(r"\d*", "", pn[0])
        pos = pn[1:4].astype(np.float64)
        struct += ase.Atom(at, position=pos)

    struct.set_cell(box + angles)
    struct.set_pbc(pbc)
    ase.io.write(f"{name}.{ft}", struct)
    print(f"{name}.{ft} - Saved! ")


def make_MOF(template, n_node_type=2, n_max_atoms=200, nodes_path="./nodes", edges_path="./edges", outdir="./outputs"):
    """
    Generate a MOF.

    Function derived from the run_template function, this function is the central use of ToBaCco.

    Parameters:
    -----------
    template : str
        Topology path.

    nodes_path : str
        Nodes path.

    edges_path : str
        Edges path.

    n_node_type : int
        Maximum number of node types (nodes with different geometries).

    n_max_atoms : int
        Maximum number of atoms allowed per structure.
    """
    print()
    print('======================================================================================')
    print('template:', template)
    print('======================================================================================')
    print()

    print("Number of node types allowed:", n_node_type)
    print(f"Number of atoms allowed: <= {n_max_atoms}")

    cat_count = 0
    for net in ct2g(template):
        cat_count += 1
        TG, start, unit_cell, TVT, TET, TNAME, a, b, c, ang_alpha, ang_beta, ang_gamma, max_le, catenation = net
        # ---------------------------------------------
        # Information:
        #   TG: Main network.
        #   start: Coordinates of the first node.
        #   unit_cell: Unit cell matrix.
        #   TVT: Set composed of tuples with the number of connections and a label
        # (example: {(4, 'Ti'), (3, 'Er'), (6, 'V')}).
        #   TET: Set composed of tuples with pairs of nodes
        # (examples: TET: {('Er', 'V'), ('Er', 'Ti')}).
        #   TNAME: template name.
        #   a, b, c, ang_alpha, ang_beta, ang_gamma: cell parameters
        #   max_le: Max lenght (float)
        #   catenation: bool
        # ---------------------------------------------
        if len(TVT) > n_node_type:
            print("Topology with number of node types greater than the one defined (%s)" % n_node_type)
            return None

        TVT = sorted(TVT, key=lambda x: x[0], reverse=True)
        print("TVT:", TVT)
        TET = sorted(TET, reverse=True)
        print("TET:", TET)
        # Node database files, with their number of connections.
        # example: [(2, '2X_Sc_D*h.cif'), (4, '4X_Sc_D4h.cif'), ...]
        node_cns = [(cncalc(node, path=nodes_path), node) for node in os.listdir(nodes_path)]
        print("node_cns:", node_cns)
        print('Number of vertices = ', len(TG.nodes()))
        print('Number of edges = ', len(TG.edges()))
        print()

        # Number of edges present.
        edge_counts = dict((data['type'], 0) for e0, e1, data in TG.edges(data=True))
        for e0, e1, data in TG.edges(data=True):
            edge_counts[data['type']] += 1

        print("edge_counts:", edge_counts)
        # if it's empty, not nodes weren't assign
        vas = vertex_assign(TG, TVT, node_cns, unit_cell, USER_SPECIFIED_NODE_ASSIGNMENT, SYMMETRY_TOL, ALL_NODE_COMBINATIONS)
        CB, CO = cycle_cocyle(TG)
        for va in vas:
            if len(va) == 0:
                print('At least one vertex does not have a building block with the correct number of connection sites.')
                print('Moving to the next template...')
                print()
                continue

        if len(CB) != (len(TG.edges()) - len(TG.nodes()) + 1):
            print('The cycle basis is incorrect.')
            print('The number of cycles in the cycle basis does not equal the rank of the cycle space.')
            print('Moving to the next tempate...')
            continue

        num_edges = len(TG.edges())
        print("num_edges:", num_edges)
        Bstar, alpha = Bstar_alpha(CB, CO, TG, num_edges)

        num_vertices = len(TG.nodes())
        print("num_vertices:", num_vertices)

        # Help to combinate edges in the assignment.
        if COMBINATORIAL_EDGE_ASSIGNMENT:
            eas = list(itertools.product([e for e in os.listdir(edges_path)], repeat=len(TET)))
        else:
            edge_files = sorted([e for e in os.listdir(edges_path)])
            eas = []
            i = 0
            while len(eas) < len(TET):
                eas.append(edge_files[i])
                i += 1
                if i == len(edge_files):
                    i = 0
            eas = [eas]
        print("eas:", eas)

        g = 0
        # Here we will test different combinations of nodes and edges
        for va in vas:
            # Selection of vertex(node) combination
            # example: [('V1', '3X_Co.cif'), ('V2', '3X_Co.cif')]
            node_elems = [bbelems(i[1], nodes_path) for i in va]
            # print("node_elems:", node_elems)
            metals = [[i for i in j if i in metal_elements] for j in node_elems]
            metals = list(set([i for j in metals for i in j]))
            # print("metals:", metals)

            v_set = [('v' + str(vname_dict[re.sub('[0-9]', '', i[0])]), i[1]) for i in va]
            v_set = sorted(list(set(v_set)), key=lambda x: x[0])
            v_set = [v[0] + '-' + v[1] for v in v_set]

            print('*****************************************************************')
            print('vertex assignment : ', v_set)
            print('*****************************************************************')
            print()

            if SINGLE_METAL_MOFS_ONLY and len(metals) != 1:
                print(v_set, 'contains no metals or multiple metal elements, no cif will be written')
                print()
                continue

            if MOFS_ONLY and len(metals) < 1:
                print(v_set, 'contains no metals, no cif will be written')
                print()
                continue

            for v in va:
                for n in TG.nodes(data=True):
                    if v[0] == n[0]:
                        n[1]['cifname'] = v[1]

            for ea in eas:
                g += 1
                print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                print('edge assignment : ', ea)
                print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                print()

                type_assign = dict((k, []) for k in sorted(TET, reverse=True))
                for k, m in zip(TET, ea):
                    type_assign[k] = m

                print("type_assign:", type_assign)
                for e in TG.edges(data=True):
                    ty = e[2]['type']
                    for k in type_assign:
                        if ty == k or (ty[1], ty[0]) == k:
                            e[2]['cifname'] = type_assign[k]

                num_possible_XX_bonds = 0
                for edge_type, cifname in zip(TET, ea):
                    if cifname == 'ntn_edge.cif':
                        factor = 1
                    else:
                        factor = 2
                    edge_type_count = edge_counts[edge_type]
                    print("edge_type_count:", edge_type_count)
                    num_possible_XX_bonds += factor * edge_type_count
                print("num_possible_XX_bonds:", num_possible_XX_bonds)
                # Here, it's where tobacco place the sbus.
                ea_dict = assign_node_vecs2edges(TG, unit_cell, SYMMETRY_TOL, template)
                all_SBU_coords = SBU_coords(TG, ea_dict, CONNECTION_SITE_BOND_LENGTH, edges_path)

                # REVISAR??
                sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma, sc_covar, Bstar_inv, max_length, callbackresults, ncra, ncca, scaling_data = scale(all_SBU_coords,a,b,c,ang_alpha,ang_beta,ang_gamma,max_le,num_vertices,Bstar,alpha,num_edges,FIX_UC,SCALING_ITERATIONS,PRE_SCALE,MIN_CELL_LENGTH,OPT_METHOD)

                print('*******************************************')
                print('The scaled unit cell parameters are : ')
                print('*******************************************')
                print('a    :', np.round(sc_a, 5))
                print('b    :', np.round(sc_b, 5))
                print('c    :', np.round(sc_c, 5))
                print('alpha:', np.round(sc_alpha, 5))
                print('beta :', np.round(sc_beta, 5))
                print('gamma:', np.round(sc_gamma, 5))
                print()

                for sc, name in zip((sc_a, sc_b, sc_c), ('a', 'b', 'c')):
                    cflag = False
                    if sc == MIN_CELL_LENGTH:
                        print('unit cell parameter', name, 'may have collapsed during scaling!')
                        print('try re-running with', name, 'fixed or a larger MIN_CELL_LENGTH')
                        print('no cif will be written')
                        cflag = True

                if cflag:
                    continue

                scaled_params = [sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma]
                sc_Alpha = np.r_[alpha[0:num_edges-num_vertices+1, :], sc_covar]
                sc_omega_plus = np.dot(Bstar_inv, sc_Alpha)

                ax = sc_a
                ay = 0.0
                az = 0.0
                bx = sc_b * np.cos(sc_gamma * pi/180.0)
                by = sc_b * np.sin(sc_gamma * pi/180.0)
                bz = 0.0
                cx = sc_c * np.cos(sc_beta * pi/180.0)
                cy = (sc_c * sc_b * np.cos(sc_alpha * pi/180.0) - bx * cx) / by
                cz = (sc_c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
                sc_unit_cell = np.asarray([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T

                scaled_coords = omega2coords(start, TG, sc_omega_plus, (sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma), num_vertices, template, g, WRITE_CHECK_FILES)
                # Here scaled node and edge place X to direction for topologies
                nvecs, evecs = scaled_node_and_edge_vectors(scaled_coords, sc_omega_plus, sc_unit_cell, ea_dict)
                # Coordinates and bond list
                try:
                    placed_nodes, node_bonds = place_nodes(
                        nvecs, CHARGES, ORIENTATION_DEPENDENT_NODES, nodes_path
                    )
                except np.linalg.LinAlgError:
                    print("convergence not achieved!")
                    return None
                # Coordinates of edges
                # Center and X--X positions
                # print("evecs:\n")
                # for ev in evecs:
                #     print(ev)
                placed_edges, edge_bonds = place_edges(evecs, CHARGES, len(placed_nodes), edges_path)
                # if RECORD_CALLBACK:
                #     vnames = '_'.join([v.split('.')[0] for v in v_set])
                #     if len(ea) <= 5:
                #         enames = '_'.join([e[0:-4] for e in ea])
                #     else:
                #         enames = str(len(ea)) + '_edges'
                #     prefix = template[0:-4] + '_' +  vnames + '_' + enames
                #     frames = scaling_callback_animation(callbackresults, alpha, Bstar_inv, ncra, ncca, num_vertices, num_edges, TG, template, g, False)
                #     write_scaling_callback_animation(frames, prefix)
                #     animate_objective_minimization(callbackresults, prefix)

                if PLACE_EDGES_BETWEEN_CONNECTION_POINTS:
                    placed_edges = adjust_edges(placed_edges, placed_nodes, sc_unit_cell)

                placed_nodes = np.c_[placed_nodes, np.array(['node' for i in range(len(placed_nodes))])]
                #######REMOVE
                # save_state(placed_nodes, [sc_a, sc_b, sc_c], [sc_alpha, sc_beta, sc_gamma], name="test_placed_nodes")
                #######

                placed_edges = np.c_[placed_edges, np.array(['edge' for i in range(len(placed_edges))])]
                #######REMOVE
                # save_state(placed_edges, [sc_a, sc_b, sc_c], [sc_alpha, sc_beta, sc_gamma], name="test_placed_edges")
                #######

                placed_all = list(placed_nodes) + list(placed_edges)
                #######REMOVE
                # save_state(placed_all, [sc_a, sc_b, sc_c], [sc_alpha, sc_beta, sc_gamma], name="test_placed_all")
                #######
                bonds_all = node_bonds + edge_bonds

                # if WRITE_CHECK_FILES:
                #     write_check_cif(template, placed_nodes, placed_edges, g, scaled_params, sc_unit_cell)

                if REMOVE_DUMMY_ATOMS:
                    # REVISAR
                    placed_all, bonds_all, nconnections = remove_Fr(placed_all, bonds_all)
                #######REMOVE
                # save_state(placed_all, [sc_a, sc_b, sc_c], [sc_alpha, sc_beta, sc_gamma], name="test_placed_all_noDummy")
                #######

                vnames = '_'.join([v.split('.')[0] for v in v_set])
                enames_list = [e[0:-4] for e in ea]
                enames_grouped = [list(edge_gr) for ind, edge_gr in itertools.groupby(enames_list)]
                enames_grouped = [(len(edge_gr), list(set(edge_gr))) for edge_gr in enames_grouped]
                enames_flat = [str(L) + '-' + '_'.join(names) for L, names in enames_grouped]
                enames = '_'.join(enames_flat)

                template_name = os.path.basename(template).split(".")[0]
                bond_check_code = ""
                if catenation:
                    cifname = template_name + '_' + vnames + '_' + enames + bond_check_code + '_' + 'CAT' + str(cat_count)  # + '.cif'
                else:
                    cifname = template_name + '_' + vnames + '_' + enames + bond_check_code  # + '.cif'

                print("Number of total atoms:", len(placed_nodes))
                if len(placed_nodes) > n_max_atoms:
                    print("Maximum number of atoms reached.")
                    return None

                if WRITE_CIF:
                    print('writing cif...')
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                    save_state(
                        placed_all,
                        [sc_a, sc_b, sc_c],
                        [sc_alpha, sc_beta, sc_gamma],
                        name=f"{outdir}/{cifname}"
                    )


                ########## For me it's not necessary.
                '''

                print(nconnections)


                print('computing X-X bonds...')
                print()
                print('*******************************************')
                print('Bond formation : ')
                print('*******************************************')

                try:
                    # REVISAR
                    fixed_bonds, nbcount, bond_check_passed = bond_connected_components(placed_all, bonds_all, sc_unit_cell, max_length, BOND_TOL, nconnections, num_possible_XX_bonds)
                except ValueError:
                    print("bond counts do not match, there is a problem with the connection site bonding algorithm")
                    return None

                break
                print('there were ', nbcount, ' X-X bonds formed')

                if bond_check_passed:
                    print('bond check passed')
                    bond_check_code = ''
                else:
                    print('bond check failed, attempting distance search bonding...')
                    fixed_bonds, nbcount = distance_search_bond(placed_all, bonds_all, sc_unit_cell, 2.5)
                    bond_check_code = '_BOND_CHECK_FAILED'
                    print('there were', nbcount, 'X-X bonds formed')
                print()

                if CHARGES:
                    fc_placed_all, netcharge, onetcharge, rcb = fix_charges(placed_all)
                else:
                    fc_placed_all = placed_all

                fixed_bonds = fix_bond_sym(fixed_bonds, placed_all, sc_unit_cell)

                if CHARGES:
                    print('*******************************************')
                    print('Charge information :                       ')
                    print('*******************************************')
                    print('old net charge                  :', np.round(onetcharge, 5))
                    print('rescaling magnitude             :', np.round(rcb, 5))

                    remove_net = choice(range(len(fc_placed_all)))
                    fc_placed_all[remove_net][4] -= np.round(netcharge, 4)
                    print('new net charge (after rescaling):', np.sum([li[4] for li in fc_placed_all]))
                    print()

                vnames = '_'.join([v.split('.')[0] for v in v_set])
                enames_list = [e[0:-4] for e in ea]
                enames_grouped = [list(edge_gr) for ind, edge_gr in itertools.groupby(enames_list)]
                enames_grouped = [(len(edge_gr), list(set(edge_gr))) for edge_gr in enames_grouped]
                enames_flat = [str(L) + '-' + '_'.join(names) for L,names in enames_grouped]
                enames = '_'.join(enames_flat)

                template_name = os.path.basename(template).split(".")[0]
                if catenation:
                    cifname = template_name + '_' + vnames + '_' + enames + bond_check_code + '_' + 'CAT' + str(cat_count) + '.cif'
                else:
                    cifname = template_name + '_' + vnames + '_' + enames + bond_check_code + '.cif'

                    print("Number of total atoms:", len(fc_placed_all))
                    if len(fc_placed_all) > n_max_atoms:
                        print("Maximum number of atoms reached.")
                        return None

                if WRITE_CIF:
                    print('writing cif...')
                    print()
                    if len(cifname) > 255:
                        cifname = cifname[0:241]+'_truncated.cif'
                    write_cif(fc_placed_all, fixed_bonds, scaled_params, sc_unit_cell, cifname, CHARGES)
                '''


def run_tobacco_serial(templates, **kwargs):
    """Run ToBaCco in serial."""
    # i = 0
    for name in templates:
        make_MOF(templates[name], **kwargs)
        # i += 1
        # if i == 2:
        #     break


def run_tobacco_parallel(templates, n_node_type=10, n_max_atoms=100):
    """Run ToBaCco in parallel."""
    poolSize = multiprocessing.cpu_count()
    print('Running parallel on', poolSize, 'processors...')
    args = [(templates[name], n_node_type, n_max_atoms) for name in templates]
    #TODO
    # def wrapper_make_MOF(template, n_node_type=n_node_type, n_max_atoms=n_max_atoms):
    #     return make_MOF(template, n_node_type=n_node_type, n_max_atoms=n_max_atoms)

    # make_MOF_partial = partial(wrapper_make_MOF, n_node_type=n_node_type)

    with Pool(processes=poolSize) as pool:
        for _ in pool.starmap(make_MOF, args):
            pass


def read_input_gaus(file):
    """Read a Gaussian .com file."""
    atoms = re.compile(r"""
        ^\s*
        (?P<atsb>\w+)\s+             # Atom name.
        (?P<x>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Z.
    """, re.X)

    xyz = []
    with open(file, "r") as XYZ:
        for line in XYZ:
            if atoms.match(line):
                m = atoms.match(line)
                xyz.append(m.groupdict())
    coord = pd.DataFrame(xyz)
    coord = coord.astype({
        "x": np.float64,
        "y": np.float64,
        "z": np.float64
    })

    return coord


def gen_geometries_metal(metal, d=0.8):
    """
    Generate geometries using a metal center.

    Parameters:
    -----------
    metal : str
        Indicates the atomic element to be used.
    """
    for pointgroup in skeleton_X:
        print(pointgroup, metal)
        sbu = build_sbu_metal_center(metal, pointgroup, d)
        parameters = extract_info(sbu)
        nX = sum(sbu.symbols == "X")
        # write cif
        write_cif_SBU(f"./nodes/{nX}X_{metal}_{pointgroup}.cif", *parameters)


def clean_folder(folder):
    """Clean folder."""
    if not os.path.exists(folder):
        os.mkdir(folder)
        print("Directory '%s' created" % folder)
    else:
        if os.path.isdir(folder):
            shutil.rmtree(folder)
            os.mkdir(folder)
            print("Directory '%s' exists and it's clean" % folder)


dftb_in = """
Geometry = vaspFormat {
<<< "MOF_PATH"
}

Driver = GeometryOptimization {
    Optimiser = Rational {}
    LatticeOpt = No
    MaxSteps = 1000
    AppendGeometries = Yes
    Convergence = {GradElem = 1E-5}
    MovedAtoms = 1:-1
    # FixAngles = No
    # FixLengths = {
    #    No No No
    # }
    # Isotropic = No
}

Hamiltonian = xTB {
  Method = "GFN2-xTB"
  kPointsAndWeights = SuperCellFolding {
    4   0   0
    0   4   0
    0   0   4
    0.5 0.5 0.5
  }
  SCC = Yes
  SCCTolerance = 1e-5
  MaxSCCIterations = 1000
  #SlaterKosterFiles = Type2FileNames {  # File names with two atom type names
  #   Prefix = "/home/ovillegas/.applications/Slater-Koster_files/3ob-3-1/"    # Prefix before first type name
  #   Separator = "-"                     # Dash between type names
  #   Suffix = ".skf"                     # Suffix after second type name
  #}
  # MaxAngularMomentum = {
  #   C = "p"
  #   O = "p"
  #   Ca = "d"
  # }
}

"""


def setting_dftb_inputs(out="./method_II_dftbplus"):
    """Configuring input files for calculations using DFTB+."""
    # Clean
    clean_folder(out)
    # Record file
    # id_mof, name, path
    with open(f"{out}/data.csv", "w") as File:
        File.write("id;name;path\n")

    for file in glob.glob("output_vasps/*.vasp"):
        name, _ = os.path.splitext(os.path.basename(file))
        idhash = name
        id_mof = hashlib.shake_256(idhash.encode("UTF-8")).hexdigest(8)
        line = "%s;%s;%s\n" % (id_mof, name, os.path.abspath(file))
        with open(f"{out}/data.csv", "a") as File:
            File.write(line)

        out_mof = f"{out}/{id_mof}"
        os.mkdir(out_mof)
        dftb_input = dftb_in.replace("MOF_PATH", os.path.abspath(file))
        with open(f"{out}/{id_mof}/dftb_in.hsd", "w") as F:
            F.write(dftb_input)
