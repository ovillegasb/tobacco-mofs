import unittest
import os
import shutil
import tobacco.tools as tools
from tobacco.ciftemplate2graph import CrystalGraph


class SettingsGENTest(unittest.TestCase):

    def setUp(self):
        topols_dict = tools.load_database()
        self.cg = CrystalGraph(topols_dict["sql"])

    def test_load_cif_reads_lines_correctly(self):
        self.assertTrue(len(self.cg.template) > 0)
        self.assertIsInstance(self.cg.template[0], str)

    def test_build_graph_creates_nodes_and_edges(self):
        G = self.cg.graph
        self.assertTrue(G.number_of_nodes() > 0)
        self.assertTrue(G.number_of_edges() > 0)


class MakeMOF2DTest(unittest.TestCase):

    def setUp(self):
        topols_dict = tools.load_database()
        self.options = {
            "template": topols_dict["sql"]
        }

    def tearDown(self):
        if os.path.exists("./nodes"):
            shutil.rmtree("./nodes")

        if os.path.exists("./edges"):
            shutil.rmtree("./edges")

    def test_make_mof_2D(self):
        metal_options = {
            "metal": "Sc",
            "pointgroup": "D4h",
            "distance": 0.8,
            "output": None
        }

        tools.gen_sbu_metal_center(**metal_options)

        edge_options = {
            "file": "file_N3.com",
            "ligand": "file_N3.com",
            "ndx_X": [0, 2],
            "output": None
        }
        tools.gen_sbu_edge(**edge_options)

        result = tools.make_MOF(**self.options)
        self.assertEqual(result, "Done!")
