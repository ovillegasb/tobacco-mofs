import numpy as np
from numpy.linalg import norm
from tobacco.bbcif_properties import bb2array, X_vecs, bbbonds, bbcharges
from tobacco.Bio import SVDSuperimposer
from tobacco import configuration
from scipy.spatial.transform import Rotation

import itertools
import re
import os
import random

# Database path
database = os.path.dirname(__file__)
nodes_path = os.path.join(database, "nodes")
edges_path = os.path.join(database, "edges")

EDGE_ROTATION = configuration.EDGE_ROTATION


def rotate_axis_XX(coord, names):
    """Rotating an edge using its X--X axis."""
    # print("Rotating the edge ...")
    # where are the dummies
    ndx_X = [i for i, name in enumerate(names) if "X" in name]
    # 2X coordinates
    ax1, ax2 = coord[ndx_X]
    # center of coordinates
    com = np.mean(coord, axis=0)
    coord -= com

    # rotation venctor
    rot_vec = ax2 - ax1
    rot_vec = rot_vec / np.linalg.norm(rot_vec)

    angles = np.arange(15, 95, 5)
    ang_deg = random.choice([-1, 1]) * random.choice(angles)
    rotacion = Rotation.from_rotvec(np.deg2rad(ang_deg) * rot_vec)
    news_coord = rotacion.apply(coord)
    news_coord += com

    return news_coord


def match_vectors(a1,a2,num):

	dist1 = [(np.linalg.norm(a1[0]-a1[i]),i) for i in range(len(a1))]
	dist2 = [(np.linalg.norm(a2[0]-a2[i]),i) for i in range(len(a2))]

	dist1.sort(key=lambda x: x[0])
	dist2.sort(key=lambda x: x[0])

	vecs1 = np.array([a1[i] for i in [dist1[j][1] for j in range(num)]])
	vecs2 = np.array([a2[i] for i in [dist2[j][1] for j in range(num)]])
	
	return vecs1,vecs2

def mag_superimpose(a1,a2):
	
	sup = SVDSuperimposer()

	a1 = np.asarray(a1)
	a2 = np.asarray(a2)
	mags = [norm(v) for v in a2]

	if len(a1) <= 7:

		min_dist = (1.0E6, 'foo', 'bar')
		
		for l in itertools.permutations(a1):

			p = np.array([m*v/norm(v) for m,v in zip(mags,l)])
			sup.set(a2,p)
			sup.run()
			rot,tran = sup.get_rotran()
			rms = sup.get_rms()

			if rms < min_dist[0]:
				min_dist = (rms,rot,tran)
	
	else:

		a1,a2 = match_vectors(a1,a2,6)
		mags = [norm(v) for v in a2]
		
		min_dist = (1.0E6, 'foo', 'bar')
		
		for l in itertools.permutations(a1):

			p = np.array([m*v/norm(v) for m,v in zip(mags,l)])
			sup.set(a2,p)
			sup.run()
			rot,tran = sup.get_rotran()
			rms = sup.get_rms()
		
			if rms < min_dist[0]:
				min_dist = (rms,rot,tran)

	return min_dist

def superimpose(a1,a2):
	
	sup = SVDSuperimposer()

	a1 = np.asarray(a1)
	a2 = np.asarray(a2)

	if len(a1) <= 7:

		min_dist = (1.0E6, 'foo', 'bar')
		
		for l in itertools.permutations(a1):

			p = np.asarray(l)
			sup.set(a2,p)
			sup.run()
			rot,tran = sup.get_rotran()
			rms = sup.get_rms()

			if rms < min_dist[0]:
				min_dist = (rms,rot,tran)
	
	else:

		a1,a2 = match_vectors(a1,a2,6)		
		min_dist = (1.0E6, 'foo', 'bar')
		
		for l in itertools.permutations(a1):

			p = np.asarray(l)
			sup.set(a2,p)
			sup.run()
			rot,tran = sup.get_rotran()
			rms = sup.get_rms()
		
			if rms < min_dist[0]:
				min_dist = (rms,rot,tran)

	return min_dist


def scaled_node_and_edge_vectors(sc_coords, sc_omega_plus, sc_unit_cell, ea_dict):
	"""Scale and place the nodes and edges."""
	# print("Vamos a desplazar y alinear los EDGES")
	# sc_coords
	# contine las posiciones de los nodos extraidos de la topologia.
	# Su numero es igual al numero de nodos en el sitema.
	nvecs = []
	evecs = []
	already_placed_edges = []
	nvecs_append = nvecs.append
	evecs_append = evecs.append
	already_placed_edges_append = already_placed_edges.append
	# print("sc_unit_cell:\n", sc_unit_cell)
	# print("sc_coords:\n", sc_coords)
	# print("ea_dict:\n", ea_dict)

	for n in sc_coords:
		# node label, cif file node, node pos, edges list []
		vertex, vcif, vfvec, indicent_edges = n
		# pos in cell (1x3)
		vcvec = np.dot(sc_unit_cell, vfvec)
		# print("Node:", vertex, vcvec, vcif)

		#
		ie = []
		ie_append = ie.append

		# print("Edges:")
		for e in indicent_edges:
			# e is a edge tuple class
			# ex: (1, ('V1', 'V2'), '2B_dmet_bzn.cif')
			# print(e)
			ind = e[0]
			positive_direction = e[1]
			ecif = e[2]

			if vertex == positive_direction[0]:
				direction = 1
				on = positive_direction[1]
			else:
				direction = -1
				on = positive_direction[0]

			dxn = ea_dict[vertex][ind][1]
			dxon = ea_dict[on][ind][1]

			ie_append((ind, direction, ecif, dxn, dxon))

		efvec = []
		ecvec = []
		efvec_append = efvec.append
		ecvec_append = ecvec.append

		# edges and its directions
		for e in ie:
			# edge by vertex
			# ex: (label_edge, label_node, edge_cif, 0.9999999996888641, 1.00000000034225)
			ind, d, ecif, dxn, dxon = e
			cs = np.dot(sc_unit_cell, sc_omega_plus[ind - 1])

			fvec = vfvec + d * sc_omega_plus[ind - 1]
			cvec = vcvec + d * cs

			ec1 = vcvec + d * dxn * (cs/np.linalg.norm(cs))
			ec2 = cvec - d * dxon * (cs/np.linalg.norm(cs))
			ecoords = np.average([ec1, ec2], axis=0)
			
			ecvec_append(cvec)
			efvec_append(fvec)

			if ind not in already_placed_edges:
				evecs_append((ind, ecif, ecoords, np.array([vcvec, cvec])))
				already_placed_edges_append(ind)

		nvecs_append((vertex, vcvec, vcif, np.asarray(ecvec)))

	return nvecs, evecs


def place_nodes(nvecs, charges, ORIENTATION_DEPENDENT_NODES, nodes_path):

	placed_nbb_coords = []
	placed_nbb_coords_extend = placed_nbb_coords.extend
	all_bonds = []
	all_bonds_extend = all_bonds.extend
	ind_seg = 0
	bbind = 1

	for n in nvecs:

		bbind = bbind + 1
		name,cvec,cif,nvec = n
		ll = 0

		for v in nvec:
			mag = np.linalg.norm(v - np.average(nvec, axis = 0))
			if mag > ll:
				ll = mag

		bbxvec = np.array(X_vecs(cif,'nodes',False))

		#if ORIENTATION_DEPENDENT_NODES:
		nbbxvec = bbxvec
		#else:
		#	nbbxvec = np.array([ll*(v / np.linalg.norm(v)) for v in bbxvec])

		min_dist,rot,tran = superimpose(nbbxvec,nvec)

		all_bb = bb2array(cif, nodes_path)
		all_coords = np.array([v[1] for v in all_bb])
		all_inds = np.array([v[0] for v in all_bb])
		chg, elem = bbcharges(cif, nodes_path)
		all_names = [o + re.sub('[A-Za-z]','',p) for o,p in zip(elem,all_inds)]

		all_names_indices = np.array([int(re.sub('[A-Za-z]','',e)) for e in all_names]) + ind_seg

		elem_dict = dict((k,'') for k in all_inds)
		for i,j in zip(all_inds, elem):
			elem_dict[i] = j

		ind_dict = dict((k,'') for k in all_inds)
		for i,j in zip(all_inds, all_names_indices):
			ind_dict[i] = j

		bonds = bbbonds(cif, nodes_path)

		anf = [str(elem_dict[n]) + str(ind_dict[n]) for n in all_inds]

		abf = []
		for b in bonds:
			b1 = str(elem_dict[b[0]]) + str(ind_dict[b[0]])
			b2 = str(elem_dict[b[1]]) + str(ind_dict[b[1]])
			abf.append([b1,b2] + b[2:])

		aff_all = np.dot(all_coords,rot) + cvec
		
		laff_all = np.c_[anf, aff_all, chg, all_inds, [bbind] * len(anf)]

		placed_nbb_coords_extend(laff_all)
		all_bonds_extend(abf)
		ind_seg = ind_seg + len(all_names)

	return placed_nbb_coords, all_bonds


def place_edges(evecs, charges, nnodes, edges_path):
	"""
	Put EDGES in place.
	
	This function applies the necessary operations to move and rotate the EDGES to the required
	position.
	"""

	placed_ebb_coords = []
	placed_ebb_coords_extend = placed_ebb_coords.extend
	all_bonds = []
	all_bonds_extend = all_bonds.extend
	ind_seg = nnodes
	bbind = -1
	print("Start to align the edges")
	print("evecs:", len(evecs))
	# This is a list of tuples like:
	# [(label, center )]
	for e in evecs:
		# e: index, cif, ecoords, evec
		# ex: (
		#         1,
		#         '2B_dmet_bzn.cif',
		#         array([3.11941231e+00, 5.40442166e+00, 1.13841349e-03]),
		#         array([[-6.23986164e-04,  7.20577549e+00,  0.00000000e+00],
		#                [ 6.23944861e+00,  3.60306783e+00,  2.27682698e-03]])
		# )
		bbind = bbind - 1
		index, cif, ecoords, evec = e
		ll = 0

		for v in evec:
			mag = np.linalg.norm(v - np.average(evec, axis=0))
			if mag > ll:
				ll = mag

		# extract X vector
		bbxvec = np.array(X_vecs(cif, edges_path, False))
		nbbxvec = np.array([ll*(v / np.linalg.norm(v)) for v in bbxvec])

		min_dist, rot, tran = superimpose(nbbxvec, evec)
		# print("min_dist:", min_dist)

		all_bb = bb2array(cif, edges_path)
		all_coords = np.array([v[1] for v in all_bb])
		all_inds = np.array([v[0] for v in all_bb])
		chg, elem = bbcharges(cif, edges_path)
		all_names = [o + re.sub('[A-Za-z]', '', p) for o, p in zip(elem, all_inds)]

		all_names_indices = np.array([int(re.sub('[A-Za-z]', '', e)) for e in all_names]) + ind_seg
		# news indices

		elem_dict = dict((k, '') for k in all_inds)
		for i, j in zip(all_inds, elem):
			elem_dict[i] = j

		ind_dict = dict((k, '') for k in all_inds)
		for i, j in zip(all_inds, all_names_indices):
			ind_dict[i] = j

		bonds = bbbonds(cif, edges_path)
		anf = [str(elem_dict[n]) + str(ind_dict[n]) for n in all_inds]

		abf = []
		for b in bonds:
			b1 = str(elem_dict[b[0]]) + str(ind_dict[b[0]])
			b2 = str(elem_dict[b[1]]) + str(ind_dict[b[1]])
			abf.append([b1, b2] + b[2:])

		# news coordinates
		aff_all = np.dot(all_coords, rot) + ecoords
		if EDGE_ROTATION:
			aff_all = rotate_axis_XX(aff_all, all_inds)

		laff_all = np.c_[anf, aff_all, chg, all_inds, [bbind] * len(anf)]

		placed_ebb_coords_extend(laff_all)
		all_bonds_extend(abf)
		ind_seg = ind_seg + len(all_names)

	return placed_ebb_coords, all_bonds
