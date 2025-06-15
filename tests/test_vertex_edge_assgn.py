import unittest
import os
import shutil
from tobacco import tools
from tobacco.bbcif_properties import cncalc
from tobacco.ciftemplate2graph import CrystalGraph
from tobacco.vertex_edge_assign import vertex_assign
from tobacco import configuration as config


class edgeAssignTest(unittest.TestCase):

    def setUp(self):
        topols_dict = tools.load_database()
        template = topols_dict["sql"]
        nodes_path = "./nodes"
        options = {
            "metal": "Sc",
            "pointgroup": "D4h",
            "distance": 0.8,
            "output": None
        }
        # Generation metal
        tools.gen_sbu_metal_center(**options)
        cg = CrystalGraph(template)
        first_net = next(iter(cg))
        first_net.TVT = sorted(first_net.TVT, key=lambda x: x[0], reverse=True)
        first_net.TET = sorted(first_net.TET, reverse=True)
        node_cns = [(cncalc(node, path=nodes_path), node) for node in os.listdir(nodes_path)]

        # Number of edges present.
        edge_counts = dict(
            (data['type'], 0) for e0, e1, data in first_net.TG.edges(data=True)
        )
        for e0, e1, data in first_net.TG.edges(data=True):
            edge_counts[data['type']] += 1

        self.options = (
            first_net.TG, first_net.TVT, node_cns, first_net.unit_cell,
            config.USER_SPECIFIED_NODE_ASSIGNMENT,
            config.SYMMETRY_TOL,
            config.ALL_NODE_COMBINATIONS
        )

    def tearDown(self):
        if os.path.exists("./nodes"):
            shutil.rmtree("./nodes")

    def test_vertex_assign(self):
        result = vertex_assign(*self.options)

        # 1. Verifica que el resultado es una lista
        self.assertIsInstance(result, list)

        # 2. Verifica que haya al menos una asignación
        self.assertGreater(len(result), 0)

        # 2. Verifica que el numero de asignacion sea 4 para sql
        self.assertEqual(len(result[0]), 4)

        first_assignment = result[0]
        # 3. Verifica que cada asignación también es una lista
        self.assertIsInstance(first_assignment, list)

        # 4. Verifica que contiene tuplas (nombre del nodo, nombre del cif)
        for item in first_assignment:
            self.assertIsInstance(item, tuple)
            self.assertEqual(len(item), 2)
            node_name, cif_file = item
            self.assertTrue(isinstance(node_name, str) and node_name.startswith("V"))
            self.assertTrue(cif_file.endswith(".cif"))
