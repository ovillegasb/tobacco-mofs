
import unittest
import os
import shutil
from ase.atoms import Atoms
import tobacco.tools as tools


class SBUMetalTest(unittest.TestCase):

    def setUp(self):
        self.options = {
            "metal": "Sc",
            "pointgroup": "D4h",
            "distance": 0.8,
            "output": None
        }

    def tearDown(self):
        if os.path.exists("./nodes"):
            shutil.rmtree("./nodes")

    def test_build_sbu_metal_center_returns_valid_atoms(self):
        sbu = tools.build_sbu_metal_center(**self.options)

        # Assert it's not empty
        self.assertGreater(
            len(sbu), 0, "The generated SBU should not be empty."
        )

        # Assert correct type
        self.assertIsInstance(
            sbu, Atoms, "The returned object should be an ASE Atoms instance."
        )

    def test_extract_info_returns_valid_parms(self):
        sbu = tools.build_sbu_metal_center(**self.options)
        parameters = tools.extract_info(sbu)

        # Assert it's not empty
        self.assertGreater(
            len(parameters), 0, "Parameters should not be empty."
        )

        # Assert correct type
        self.assertIsInstance(
            parameters,
            tuple,
            "The returned object should be an tupple instance."
        )

    def test_gen_name_from_sbu_valid(self):
        sbu = tools.build_sbu_metal_center(**self.options)
        result = tools.gen_name_from_sbu(sbu, self.options["pointgroup"])
        self.assertEqual(result, "ScX4_D4h.cif")

    def test_node_directory(self):
        result = tools.node_directory()
        self.assertEqual(result, "./nodes")

    def test_node_dir_file(self):
        output = "test.cif"
        result = tools.node_directory(output=output)
        self.assertEqual(result, "./nodes/test.cif")

    def test_node_dir_path(self):
        output = "nodes/test.cif"
        result = tools.node_directory(output=output)
        self.assertEqual(result, "./nodes/test.cif")

    def test_metal_sbu_gen(self):
        result = tools.gen_sbu_metal_center(**self.options)
        # self.assertTrue(os.path.exists("./nodes/ScX4_D4h.cif"))

        self.assertEqual(result, "Done!")

    def test_gen_geometries_metal(self):
        result = tools.gen_geometries_metal(**self.options)
        self.assertEqual(result, "Done!")
