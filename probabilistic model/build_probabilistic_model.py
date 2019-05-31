
# ===========================================================================================
#
# This code generates the probabilstic model of connectivity from average properties of 
# a number m of generated anatomical connectomes. This file returns the probabilistic matrix
# of connectivity P_m.txt and the average rostro-caudal positions of neurons pos_m.txt. 
# Additionally, the code saves the adjacency matrixes of all the anatomical connectome used
# and the RC positions of the neuron is the realization in the directory 
# "anatomical adjacency matrixes". The user should specify the path of sample anatomical 
# connectomes in the variable path. If the path is not specified by the user path is 
# "../anatomical model/connectome files". 
#
# ===========================================================================================

import numpy as np
import random

n=1382 # total number of neurons
m=10 # number of connectomes 
path='../anatomical model/connectome files/' # path of the sample connetomes

P=np.zeros((n,n)) # probability matrix
pos=np.zeros((n,1)) # vector of average positions


class Cell:
	'''
	This class expresses the properties of each neuron in the model: the neuron id, type, 
	relative id, rostro-caudal position, the id of the neurons from where the neuron
	receives connections and the cell body position (1 = left, 2 = right).
	'''
	def __init__(self, id, type_id, relative_id, pos, body_side=None):
		self.id = id
		self.type_id = type_id
		self.relative_id = relative_id
		self.pos = pos
		self.incoming_connections = []
		self.body_side = body_side 

def get_cell(cell_list, id):
	'''
	check if id is cellist
	'''
	for c in cell_list:
		if c.id == id:
			return c

def load_cells(cell_file_path, connectome_file_path):
	'''
	This routine reads the DendriteGrad.txt and inc_connectGrad.txt that 
	are returned from the anatomical model of connectivity and generates
	a list of cells Cell objects with the data specified in the files
	'''
	cell_file = open(cell_file_path, "rt")
	type_counts = {}
	cell_list = []
	for l in cell_file.readlines():
		toks = l.split()
		id = int(float(toks[0])) - 1
		type_id = int(float(toks[1])) - 1
		if type_id not in type_counts:
			type_counts[type_id] = 0
		else:
			type_counts[type_id] += 1

		relative_id = type_counts[type_id]
		pos = float(toks[2])
		cell_list.append(Cell(id, type_id, relative_id, pos))
	cell_file.close()

	for c in cell_list:
		c.body_side = 1 if c.id < len(cell_list)/2 else 2

	connectome_file = open(connectome_file_path)
	connectome_file.readline() # Skip first line
	line = 1
	for l in connectome_file.readlines():
		cell = get_cell(cell_list, line-1)
		if cell is None:
			raise Exception("%s:%d: no cell with ID %d." %(connectome_file_path, line, line))

		cols = [int(s) for s in l.split()]
		for i in range((cols[0] - 1) / 2):
			inc_id = cols[1 + i*2] - 1
			inc_type = cols[2 + i*2] - 1
			inc_cell = get_cell(cell_list, inc_id)
			if inc_cell is None:
				raise Exception("%s:%d: can't find cell ID %d." % (connectome_file_path, line, inc_id+1))
			if inc_cell.type_id != inc_type:
				raise Exception("%s:%d: connectome and cell list specify different cell types (%d <> %d)." %(connectome_file_path, line, inc_cell.type_id, inc_type))
	
			cell_list[line-1].incoming_connections.append(inc_cell)
	
		line += 1
	connectome_file.close()
	
	return cell_list


# ===== BUILD THE PROBABILITY MATRIX BETWEEN EACH CELL TYPE AND EACH CELL AVERAGE RC POSITIONS =====
for connectome_id in range(1,m+1):
	print connectome_id
	# load sample connectome i
	cellist=load_cells(path+'DendriteGrad'+str(connectome_id)+'.txt',path+'inc_connectGrad'+str(connectome_id)+'.txt')
	idx=[x.type_id for x in cellist]
	idx=np.array(idx)
	
	vec_idx=[]
	for i in xrange(7):
		v=np.where(idx==i)
		vec_idx.append(v[0].tolist())
	vec_idx=np.concatenate((vec_idx)).tolist()

	# create new cellist with ordered RC positions
	new_cellist=[cellist[ind] for ind in vec_idx]
	for i in xrange(n):
		new_cellist[i].id=i
	
	A=np.zeros((n,n))
	pos_A=[]
	for i in xrange(n):
		inc_conn=[y.id for y in new_cellist[i].incoming_connections] # calculate unique vector of incoming connections (avoid multiple connections from pre-post neurons)
		inc_conn=list(set(inc_conn)) 
		A[i,inc_conn]=+1
		P[i,inc_conn]=P[i,inc_conn]+1

		# mean values update
		pos_A.append(new_cellist[i].pos)
		pos[i]+=new_cellist[i].pos
	np.save("anatomical adjacency matrixes/A_connectome"+str(connectome_id), A.transpose())
	np.save("anatomical adjacency matrixes/pos"+str(connectome_id), pos_A)

	# additionally save 
	#import scipy.io as sio
	#sio.savemat("anatomical adjacency matrixes/A_connectome"+str(connectome_id), {'A':A.transpose()})


P=P/m
pos=np.concatenate(pos/m).tolist()

np.save("anatomical adjacency matrixes/P",P)
np.save("anatomical adjacency matrixes/pos",pos)




