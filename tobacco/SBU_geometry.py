from tobacco.bbcif_properties import calc_edge_len
import numpy as np
import os

# Database path
database = os.path.dirname(__file__)


def SBU_coords(TG, ea_dict, csbl, edges_path=None):
    """
    Generate 3D connection vectors for each node (SBU) in the topology graph.

    Args:
        TG (networkx.MultiDiGraph): Topology graph where nodes are SBUs and
        edges represent connections.
        ea_dict (dict): Dictionary mapping each node and edge index to a
        vector and offset values.
                        Format: {vertex: {index: (label, dx, vector)}}
        csbl (float): Connection Site Bond Length. Acts as a buffer/extension
        distance used in the bond.
        edges_path (str): Path to the directory containing edge CIF files. If
        None, defaults to 'database/edges'.

    Returns:
        list: List of tuples (vertex, vectors), where 'vectors' is a list of
        [index, scaled_vector]
    """

    SBU_coords = []
    SBU_coords_append = SBU_coords.append

    if edges_path is None:
        edges_path = os.path.join(database, "edges")

    for node in TG.nodes(data=True):
        vertex = node[0]  # The node label (e.g., 'M1' or 'L3')
        xvecs = []
        xvecs_append = xvecs.append

        for e0, e1, edict in TG.edges(data=True):
            if vertex in (e0, e1):
                ecif = edict['cifname']             # Name of the edge file
                positive_direction = edict['pd']    # Directional (source, target)
                ind = edict['index']                # Edge index
                length = calc_edge_len(ecif, edges_path)  # Edge length from CIF

                # Determine direction of the vector based on vertex position in the edge
                if vertex == positive_direction[0]:
                    direction = 1
                    ov = positive_direction[1]  # Opposite vertex
                else:
                    direction = -1
                    ov = positive_direction[0]

                xvecname, dx_v, xvec = ea_dict[vertex][ind]
                dx_ov = ea_dict[ov][ind][1]  # Offset from opposite vertex

                # Compute total connection vector length
                if length < 0.1:  # Dummy or short bond
                    total_length = dx_v + dx_ov + csbl
                else:
                    # Actual bond + padding on both sides
                    # total_length = dx_v + dx_ov + length + 2 * csbl
                    total_length = length + 2 * csbl

                # Normalize vector and scale with total length and direction
                svec = (xvec / np.linalg.norm(xvec)) * total_length * direction
                xvecs_append([ind, svec])

        SBU_coords_append((vertex, xvecs))

    return SBU_coords
