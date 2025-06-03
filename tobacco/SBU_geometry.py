from tobacco.bbcif_properties import calc_edge_len
import numpy as np
import os

# Database path
database = os.path.dirname(__file__)

def SBU_coords(TG, ea_dict, csbl, edges_path=None):

	SBU_coords = []
	SBU_coords_append = SBU_coords.append

	if edges_path is None:
		edges_path = os.path.join(database, "edges")

	for node in TG.nodes(data=True):

		vertex = node[0]
		xvecs = []
		xvecs_append = xvecs.append

		for e0, e1, edict in TG.edges(data=True):

			if vertex in (e0,e1):

				ecif = edict['cifname']
				positive_direction = edict['pd']
				ind = edict['index']
				length = calc_edge_len(ecif, edges_path)

				if vertex == positive_direction[0]:
					direction = 1
					ov = positive_direction[1]
				else:
					direction = -1
					ov = positive_direction[0]
				
				xvecname,dx_v,xvec = ea_dict[vertex][ind]
				dx_ov = ea_dict[ov][ind][1]

				if length < 0.1:
					total_length = dx_v + dx_ov + csbl
				else:
					total_length = dx_v + dx_ov + length + 2*csbl
				
				svec = (xvec/np.linalg.norm(xvec)) * total_length * direction
				xvecs_append([ind, svec])

		SBU_coords_append((vertex, xvecs))

	return SBU_coords
