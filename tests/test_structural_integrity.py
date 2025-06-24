import unittest
from tobacco import tools
import os
import shutil
from ase.io import read
from ase.geometry.analysis import Analysis
from ase.geometry import get_distances
from ase.neighborlist import neighbor_list
from collections import defaultdict
import numpy as np


def gen_input_com():
    content = """# hf/3-21g

N3 ligand

0 2
 N                  0.94974703   -0.83333332    0.00902855
 N                 -0.32425297   -0.83333332    0.00902855
 N                 -1.59825297   -0.83333332    0.00902855

    """
    with open("N3.com", "w") as com:
        com.write(content)


def gen_input_ion_xyz():
    content = """1
Properties=species:S:1:pos:R:3 pbc="F F F"
Na       0.00000000       0.00000000       0.00000000"""
    with open("Na.xyz", "w") as com:
        com.write(content)


def gen_input_ion_NH4_xyz():
    content = """5

 H                  1.35665441   -0.80881001    0.00000000
 H                  1.35667283    0.70439819    0.87365150
 H                  1.35667283    0.70439819   -0.87365150
 H                 -0.07000001    0.20001318    0.00000000
 N                  0.99999999    0.20000000    0.00000000"""
    with open("NH4.xyz", "w") as com:
        com.write(content)


def min_distance_from_rdf(atoms, el1="Sc", el2="N", rmax=10.0, nbins=200, threshold=1e-5):
    """
    Return min value between el1 and el2 from RDF.
    """
    analysis = Analysis(atoms)

    values, distances = analysis.get_rdf(
                            rmax=rmax, nbins=nbins, elements=(el1, el2),
                            return_dists=True
                        )[0]

    idx = next(i for i, v in enumerate(values) if v > threshold)
    return distances[idx]


def get_pairwise_distances(atoms, cutoff=10.0, pairs=[('Sc', 'N'), ('Sc', 'Sc')]):
    # Obtener listas de pares de índices y sus distancias (PBC considerado)
    i_list, j_list, distances = neighbor_list('ijd', atoms, cutoff)

    # Crear índice inverso: índice -> símbolo
    symbols = atoms.get_chemical_symbols()

    # Diccionario para guardar distancias por tipo de par
    distances_by_pair = defaultdict(list)

    for i, j, d in zip(i_list, j_list, distances):
        s1, s2 = symbols[i], symbols[j]
        pair = tuple(sorted((s1, s2)))

        if pair in [tuple(sorted(p)) for p in pairs]:
            distances_by_pair[pair].append(d)

    return distances_by_pair


class GeneratedMOFIntegrityTest(unittest.TestCase):

    def setUp(self):
        self.topols_dict = tools.load_database(db_path="../porE/data_ions/")

        # self.options = {
        #     "n_max_atoms": 100,
        #     "n_node_type": 1,
        #     "ion": "./NH4.xyz",
        #     "templates_path": "",
        #     "n_ions": 1
        # }

        # metal_options = {
        #     "metal": "Sc",
        #     "pointgroup": "Oh",
        #     "distance": 0.8,
        #     "output": None
        # }

        # tools.gen_sbu_metal_center(**metal_options)
        # tools.gen_geometries_metal(**metal_options)

        # gen_input_com()
        # edge_options = {
        #     "file": "N3.com",
        #     "ligand": "N3.com",
        #     "ndx_X": [0, 2],
        #     "output": None
        # }
        # tools.gen_sbu_edge(**edge_options)

        # gen_input_ion_NH4_xyz()

    #def tearDown(self):
    #    if os.path.exists("./nodes"):
    #        shutil.rmtree("./nodes")
    #
    #    if os.path.exists("./edges"):
    #        shutil.rmtree("./edges")
    #
    #    # if os.path.exists("./outputs"):
    #    #     shutil.rmtree("./outputs")
    #
    #    if os.path.isfile("N3.com"):
    #        os.remove("N3.com")
    #
    #    if os.path.isfile("NH4.xyz"):
    #        os.remove("NH4.xyz")

    def test_inner_distances(self):
        topols_test = [
            # ("utb", "./outputs/utb_v1-ScX3_D3h_1-N3_1-NH4.cif"),
            # ("crb", "./outputs/crb_v1_ScX4_Td_1-N3_1-NH4.cif"),
            ("irl", "./outputs/irl_v1-ScX4_D4h_1-N3_1-NH4.cif"),
            ("lon", "./outputs/lon_v1-ScX4_Td_1-N3_1-NH4.cif"),
            ("nbo", "./outputs/nbo_v1-ScX4_D4h_1-N3_1-NH4.cif"),
            ("pbz", "./outputs/pbz_v1-ScX3_D3h_1-N3_1-NH4.cif"),
            ("pcb", "./outputs/pcb_v1-ScX4_Td_1-N3_1-NH4.cif"),
            # ("dft", "./outputs/dft_v1-ScX4_Td_1-N3_1-NH4.cif"),  # FAIL
            # ("atn", "./outputs/atn_v1-ScX4_Td_1-N3_1-NH4.cif"),  # FAIL
            # ("nboa", "./outputs/nboa_v1-ScX3_D3h_1-N3_1-NH4.cif"),  # FAIL
            # ("qtza", "./outputs/qtza_v1-ScX4_Td_1-N3_1-NH4.cif"),  # FAIL
            # ("ukf", "./outputs/ukf_v1-ScX4_Td_1-N3_1-NH4.cif"),  # FAIL
        ]

        for i, (topol, output) in enumerate(topols_test):
            with self.subTest(run=i+1):
                print("Topology:", topol)
                # template = self.topols_dict[topol]
                # result = tools.make_MOF(
                #     template=template,
                #     **self.options
                # )
                # self.assertEqual(result, "Done!")

                struct = read(output)

                # # Calculates minimum cell dimension
                # cell_lengths = struct.get_cell_lengths_and_angles()[:3]  # solo longitudes
                # min_length = min(cell_lengths)
                # safe_rmax = min_length / 2.01  # leaves margin
                # ScN_distance = min_distance_from_rdf(
                #     struct, el1="Sc", el2="N", rmax=safe_rmax
                # )
                # print(f"Minimum distance Sc–N ≈ {ScN_distance:.3f} Å")

                # dists, _ = get_distances(
                #     struct.get_positions(),  # posiciones cartesianas
                #     struct.get_positions(),
                #     cell=struct.get_cell(),
                #     pbc=struct.get_pbc(),
                # )
                # print(dists.shape)
                # print(dists)

                dist_dict = get_pairwise_distances(
                    struct, cutoff=10.0, pairs=[('Sc', 'N')]
                )

                # Show results
                for pair, dlist in dist_dict.items():
                    print(f"\n{pair[0]}–{pair[1]}: {len(dlist)} pairs found")
                    print(f"  Minimum distance: {min(dlist):.3f} Å")
                    print(f"  Maximum distance: {max(dlist):.3f} Å")

                self.assertGreaterEqual(np.round(min(dlist), 2), 2.00)
