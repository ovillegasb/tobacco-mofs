
import unittest
import os
import shutil
from ase.atoms import Atoms
import tobacco.tools as tools


class SBUedgesTest(unittest.TestCase):

    def setUp(self):
        self.options = {
            "file": "file_N3.com",
            "ligand": "file_N3.com",
            "ndx_X": [0, 2],
            "output": None
        }

    def tearDown(self):
        if os.path.exists("./edges"):
            shutil.rmtree("./edges")

    def _gen_input_file_com(self):
        pass

    def test_edge_directory(self):
        result = tools.edge_directory()
        self.assertEqual(result, "./edges")

    def test_edge_dir_file(self):
        output = "test.cif"
        result = tools.edge_directory(output=output)
        self.assertEqual(result, "./edges/test.cif")

    def test_edge_dir_path(self):
        output = "edges/test.cif"
        result = tools.edge_directory(output=output)
        self.assertEqual(result, "./edges/test.cif")

    def test_build_sbu_valid_from_file_com(self):
        sbu = tools.build_sbu_from_gaus(**self.options)

        # Assert it's not empty
        self.assertGreater(
            len(sbu), 0, "The generated SBU should not be empty."
        )

        # Assert correct type
        self.assertIsInstance(
            sbu, Atoms, "The returned object should be an ASE Atoms instance."
        )

    def test_extract_info_returns_valid_parms(self):
        sbu = tools.build_sbu_from_gaus(**self.options)
        ndx_X = [0, 2]
        parameters = tools.extract_info(sbu, ndx_X)

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

    def test_gen_edge_from_file_com_valid(self):
        result = tools.gen_sbu_edge(**self.options)
        self.assertEqual(result, "Done!")
