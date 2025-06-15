from __future__ import print_function
import re
import os
import numpy as np
import networkx as nx

# Database path
database = os.path.dirname(__file__)

vertices = ('V' , 'Er', 'Ti', 'Ce', 'S',
            'H' , 'He', 'Li', 'Be', 'B',
            'C' , 'N' , 'O' , 'F' , 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P',
            'Cl', 'Ar', 'K' , 'Ca', 'Sc',
            'Cr', 'Mn', 'Fe', 'Co', 'Ni')
pi = np.pi

def isfloat(value):
    """
        determines if a value is a float
    """
    try:
        float(value)
        return True
    except ValueError:
        return False

def nn(string):
    return re.sub('[^a-zA-Z]','', string)

def nl(string):
    return re.sub('[^0-9]','', string)

def isvert(line):
    """
        identifies coordinates in CIFs
    """
    if len(line) >=5: 
        if nn(line[0]) in vertices and line[1] in vertices and False not in map(isfloat,line[2:5]):
            return True
        else:
            return False
    
def isedge(line):
    """
        identifies bonding in cifs
    """
    if len(line) >=5:
        if nn(line[0]) in vertices and nn(line[1]) in vertices and isfloat(line[2]) and line[-1] in ('S', 'D', 'T', 'A'):
            return True
        else:
            return False

def PBC3DF(c1, c2):

    diffa = c1[0] - c2[0]
    diffb = c1[1] - c2[1]
    diffc = c1[2] - c2[2]

    if diffa > 0.5:
        c2[0] = c2[0] + 1.0
    elif diffa < -0.5:
        c2[0] = c2[0] - 1.0
    
    if diffb > 0.5:
        c2[1] = c2[1] + 1.0
    elif diffb < -0.5:
        c2[1] = c2[1] - 1.0
 
    if diffc > 0.5:
        c2[2] = c2[2] + 1.0
    elif diffc < -0.5:
        c2[2] = c2[2] - 1.0

    return c2


class SubgraphData:
    def __init__(self, graph, start, unit_cell, cn_types, edge_types,
                 cifname, cell_params, max_le, catenation):
        self.TG = graph
        self.start = start
        self.unit_cell = unit_cell
        self.TVT = cn_types
        self.TET = edge_types
        self.TNAME = cifname
        self.cell_params = cell_params  # dict with:
        #                               aL, bL, cL, alpha, beta, gamma
        self.max_le = max_le
        self.catenation = catenation

    def __repr__(self):
        return "<SubgraphData({}, nodes={}, edges={})>".format(
            self.cifname,
            len(self.graph.nodes),
            len(self.graph.edges))


class CrystalGraph:
    def __init__(self, path):
        self.path = path
        self.cifname = os.path.basename(path).split(".")[0]
        self.graph = nx.MultiGraph()
        self.unit_cell = None
        self.start = None
        self.max_le = 1.0e6
        self.cell_params = {}
        self.catenation = False
        self.subgraphs = []
        self.cn_types = set()
        self.edge_types = set()
        self._load_cif()
        self._build_graph()
        self._split_subgraphs()

    def _load_cif(self):
        with open(self.path, 'r') as template:
            self.template = list(filter(None, template.read().split('\n')))

    def _build_graph(self):
        G = self.graph
        edge_exist = False
        types = []
        aae = []
        for line in self.template:
            s = line.split()
            if '_cell_length_a' in line:
                self.cell_params['aL'] = float(s[1])
            if '_cell_length_b' in line:
                self.cell_params['bL'] = float(s[1])
            if '_cell_length_c' in line:
                self.cell_params['cL'] = float(s[1])
            if '_cell_angle_alpha' in line:
                self.cell_params['alpha'] = float(s[1])
            if '_cell_angle_beta' in line:
                self.cell_params['beta'] = float(s[1])
            if '_cell_angle_gamma' in line:
                self.cell_params['gamma'] = float(s[1])

        # Build unit cell
        aL = self.cell_params['aL']
        bL = self.cell_params['bL']
        cL = self.cell_params['cL']
        alpha = self.cell_params['alpha']
        beta = self.cell_params['beta']
        gamma = self.cell_params['gamma']

        ax, ay, az = aL, 0.0, 0.0
        bx = bL * np.cos(gamma * pi / 180.0)
        by = bL * np.sin(gamma * pi / 180.0)
        bz = 0.0
        cx = cL * np.cos(beta * pi / 180.0)
        cy = (cL * bL * np.cos(alpha * pi / 180.0) - bx * cx) / by
        cz = (cL**2 - cx**2 - cy**2)**0.5
        self.unit_cell = np.array([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T
        ne = 0

        for line in self.template:
            s = line.split()
            if isvert(s):
                ty = re.sub('[^a-zA-Z]', '', s[0])
                types.append(ty)
                f_nvec = np.array(list(map(float, s[2:5])))
                c_nvec = np.dot(self.unit_cell, f_nvec)
                G.add_node(
                    s[0],
                    type=ty,
                    fcoords=f_nvec,
                    ccoords=c_nvec,
                    cn=[],
                    cifname=[]
                )
            if isedge(s):
                edge_exist = True
                if '_' in s[3]:
                    lbl = np.array(list(map(int, s[3].split('_')[1]))) - 5
                elif s[3] == '.':
                    lbl = np.array([0, 0, 0])
                else:
                    raise ValueError(
                        "Unrecognized bond translational \
                            symmetries in " + self.cifname
                    )
                nlbl = -1 * lbl
                if (
                    (s[0], s[1], *lbl) not in aae and
                    (s[1], s[0], *lbl) not in aae and
                    (s[0], s[1], *nlbl) not in aae and
                    (s[1], s[0], *nlbl) not in aae
                ):
                    ne += 1
                    aae.append((s[0], s[1], *lbl))
                    v1 = G.nodes[s[0]]['fcoords']
                    v2 = G.nodes[s[1]]['fcoords'] + lbl
                    ef_coords = np.average(np.array([v1, v2]), axis=0)
                    ec_coords = np.dot(self.unit_cell, ef_coords)
                    cdist = np.linalg.norm(np.dot(self.unit_cell, v1 - v2))
                    le = float(s[2])
                    self.max_le = min(self.max_le, cdist)
                    G.add_edge(
                        s[0], s[1],
                        key=(len(aae), *lbl),
                        label=lbl,
                        length=le,
                        fcoords=ef_coords,
                        ccoords=ec_coords,
                        index=ne,
                        pd=(s[0], s[1])
                    )

        if not edge_exist:
            raise ValueError(f'No edges defined in template: {self.cifname}')

    def _split_subgraphs(self):
        S = [
            self.graph.subgraph(c).copy() for c in nx.connected_components(self.graph)
        ]
        if len(S) > 1:
            self.catenation = True

        for sub in S:
            SG = nx.MultiGraph()
            for n, data in sub.nodes(data=True):
                cn = self.graph.degree(n)
                ty = re.sub('[0-9]', '', n)
                self.cn_types.add((cn, ty))
                SG.add_node(n, **data)
                if self.start is None:
                    self.start = data['fcoords']

            for count, (e0, e1, key, data) in enumerate(
                sub.edges(keys=True, data=True), 1
            ):
                key = (count,) + key[1:]
                l = sorted([
                    re.sub('[^a-zA-Z]', '', e0), re.sub('[^a-zA-Z]', '', e1)
                ])
                self.edge_types.add(tuple(l))
                SG.add_edge(e0, e1, key=key, type=tuple(l), **data)

            self.subgraphs.append(SG)

    def __iter__(self):
        """
        Allow iteration over the generated subgraphs.
        Yields:
            SubgraphData: Encapsulated graph information.
        """
        for sg in self.subgraphs:
            yield SubgraphData(
                graph=sg,
                start=self.start,
                unit_cell=self.unit_cell,
                cn_types=self.cn_types,
                edge_types=self.edge_types,
                cifname=self.cifname,
                cell_params=self.cell_params,
                max_le=self.max_le,
                catenation=self.catenation
            )


def node_vecs(node, G, unit_cell, label):
    edge_coords = []
    edge_coords_append = edge_coords.append

    for e in G.edges(data=True):

        edict = e[2]
        positive_direction = edict['pd']
        lbl = edict['label']
        ind = edict['index']

        if node in e[0:2]:

            if node == positive_direction[0]:
                vec = edict['ccoords']
            else:
                vec = np.dot(unit_cell, -1 * lbl + edict['fcoords'])
            if label:
                edge_coords_append([ind, vec])
            else:
                edge_coords_append(vec)

    if label:
        ec_com = np.average(np.asarray([v[1] for v in edge_coords]), axis=0)
    else:
        ec_com = np.average(edge_coords, axis = 0)
    if label:
        shifted_edge_coords = [[v[0], v[1] - ec_com] for v in edge_coords]
    else:
        shifted_edge_coords = [vec - ec_com for vec in edge_coords]

    return shifted_edge_coords

def edge_vecs(edge, G, unit_cell):

    for e in G.edges(data=True):
        edict = e[2]
        if edict['index'] == edge:
            s,e = e[0:2]
            ccoords = edict['ccoords']
            v1 = G.nodes[s]['ccoords']
            v2 = G.nodes[e]['ccoords']

            return [v1 - ccoords, v2 - ccoords]
