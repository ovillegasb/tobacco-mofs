import re
import networkx as nx


def remove_Fr(placed_all, bonds_all):
	"""Remove atoms dummies using Fr as center."""
	# print(placed_all[0])
	G = nx.Graph()
	for n in placed_all:
		G.add_node(n[0])
	for l in bonds_all:
		G.add_edge(l[0],l[1])

	# print("G:", G)
	# print()
	
	nconnections = []
	for l in placed_all:
		# print(l)
	
		nbors = list(G.neighbors(l[0]))
		# print(nbors)
		
		if 'Fr' in ''.join(nbors) and 'Fr' not in l[0]:
			# print("Method 1")
			count = len([nbor for nbor in nbors if 'Fr' in nbor])
			# print(count)
			l[-3] = 'X' + re.sub('[^0-9]','',l[-3])
			# print(l)
			nconnections.append([l[0], count])

			# elif 'Fr' in ''.join(nbors) and len(nbors) == 1 and 'Fr' in l[0]:
			# print("Method 2")
			# #TODO
			# print("ELIF")
			# print(l[0], ":", ''.join(nbors))
			# exit()
			# count = 2
			# l[-3] = 'X' + re.sub('[^0-9]','',l[-3])
			# for bX in nbors:
			# 	print(list(G.neighbors(bX)))
			# # nconnections.append([l[0],count])

			# else:
			# #TODO
			# print("Method 2")
			# print(l[0], ":", ''.join(nbors))
			# exit()
			# count = 2
			# l[-3] = 'X' + re.sub('[^0-9]','',l[-3])
			# for bX in nbors:
			# 	print(list(G.neighbors(bX)))
			# # nconnections.append([l[0],count])
	# print(nconnections)

	new_placed_all = []
	new_placed_all_append = new_placed_all.append

	new_bonds_all = []
	new_bonds_all_append = new_bonds_all.append
	# print("new place")
	for l in placed_all:
		# print(l)
		if re.sub('[0-9]','',l[0]) != 'Fr':
			new_placed_all_append(l)

	# print(new_placed_all[0])

	# print("new_bonds_all")
	for l in bonds_all:
		# print(l)
		if re.sub('[0-9]','',l[0]) != 'Fr' and re.sub('[0-9]','',l[1]) != 'Fr':
			new_bonds_all_append(l)
		## elif re.sub('[0-9]','',l[0]) == 'Fr' and re.sub('[0-9]','',l[1]) == 'Fr':
		## 	new_bonds_all_append(l)

	# print(new_bonds_all)

	return (new_placed_all, new_bonds_all, nconnections)

	