import unittest
import os
import shutil
import random
from tobacco import tools


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


class MakeMOF3D_ION_Test(unittest.TestCase):
    # def setUp(self):
    #     topols_dict = tools.load_database()
    #     self.options = {
    #         "template": topols_dict["pcu"]
    #     }

    def tearDown(self):
        if os.path.exists("./nodes"):
            shutil.rmtree("./nodes")

        if os.path.exists("./edges"):
            shutil.rmtree("./edges")

        if os.path.exists("./outputs"):
            shutil.rmtree("./outputs")

        if os.path.isfile("N3.com"):
            os.remove("N3.com")

        if os.path.isfile("NH4.xyz"):
            os.remove("NH4.xyz")

        if os.path.isfile("Na.xyz"):
            os.remove("Na.xyz")

    def test_load_database_new_path(self):
        topols_dict = tools.load_database(db_path="./data")
        self.assertGreaterEqual(len(topols_dict), 1)

    def test_extract_he_coords(self):
        topols_dict = tools.load_database(db_path="./data")
        coords = tools.extract_kr_coords_from_cif(topols_dict["pcu"])
        print(coords)

    def test_make_mof_3D(self):
        metal_options = {
            "metal": "Sc",
            "pointgroup": "Oh",
            "distance": 0.8,
            "output": None
        }

        tools.gen_sbu_metal_center(**metal_options)

        gen_input_com()
        edge_options = {
            "file": "N3.com",
            "ligand": "N3.com",
            "ndx_X": [0, 2],
            "output": None
        }
        tools.gen_sbu_edge(**edge_options)

        gen_input_ion_xyz()
        topols_dict = tools.load_database(db_path="./data")
        result = tools.make_MOF(
            template=topols_dict["pcu"],
            ion="./Na.xyz",
            templates_path="",
            n_ions=8
        )

        self.assertEqual(result, "Done!")

    def test_make_mof_3D_ion_NH4(self):
        metal_options = {
            "metal": "Sc",
            "pointgroup": "Oh",
            "distance": 0.8,
            "output": None
        }

        tools.gen_sbu_metal_center(**metal_options)

        gen_input_com()
        edge_options = {
            "file": "N3.com",
            "ligand": "N3.com",
            "ndx_X": [0, 2],
            "output": None
        }
        tools.gen_sbu_edge(**edge_options)

        gen_input_ion_NH4_xyz()
        topols_dict = tools.load_database(db_path="./data")
        result = tools.make_MOF(
            template=topols_dict["pcu"],
            ion="./NH4.xyz",
            templates_path="",
            n_ions=3
        )

        self.assertEqual(result, "Done!")


class MakeMOF_ION_ParallelTest(unittest.TestCase):

    def setUp(self):
        topols_dict = tools.load_database(db_path="../porE/data_ions/")
        templates = {top: topols_dict[top] for top in random.sample(list(topols_dict.keys()), 16)}

        self.options = {
            "templates": templates,
            "n_max_atoms": 200,
            "n_node_type": 1,
            "ion": "./NH4.xyz",
            "templates_path": "",
            "n_ions": 3
        }

    def tearDown(self):
        if os.path.exists("./nodes"):
            shutil.rmtree("./nodes")

        if os.path.exists("./edges"):
            shutil.rmtree("./edges")

        if os.path.exists("./outputs"):
            shutil.rmtree("./outputs")

        if os.path.isfile("N3.com"):
            os.remove("N3.com")

        if os.path.isfile("NH4.xyz"):
            os.remove("NH4.xyz")

    def test_make_mof_3D(self):
        metal_options = {
            "metal": "Sc",
            "pointgroup": "Oh",
            "distance": 0.8,
            "output": None
        }

        # tools.gen_sbu_metal_center(**metal_options)
        tools.gen_geometries_metal(**metal_options)

        gen_input_com()
        edge_options = {
            "file": "N3.com",
            "ligand": "N3.com",
            "ndx_X": [0, 2],
            "output": None
        }
        tools.gen_sbu_edge(**edge_options)

        gen_input_ion_NH4_xyz()

        result = tools.run_tobacco_parallel(**self.options)
        self.assertEqual(result, "Done!")
