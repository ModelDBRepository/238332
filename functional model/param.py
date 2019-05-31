
'''
This file defines a SimulationParameter class containing all the parameters used in the functional model 
'''

import numpy as np
import numpy.random as rnd

class SimulationParams(object):
	def __init__(self):
		# ===== NUMBER OF CELLS AND CELL TYPES SPECIFICATION (UNIVERSAL ORDERING) =====
		self.path = "../probabilistic model/anatomical adjacency matrixes/" # specify the path where the connectome data and cell positions are saved
		# load matrix P and positions if probabilistic connectome
		self.P = np.load(self.path+'P.npy')
		self.P = self.P.transpose()
		self.pos = np.load(self.path+"pos.npy").tolist() # list of RC positions
		self.cell_types = ["rb","dlc","ain","cin","din","mn","dla"] # cell types
		self.tot_num_x_type = [0, 126, 104, 136, 384, 236, 338, 58] # total number of cells in the whole connectome 
		self.num_types = [self.tot_num_x_type[i]+sum(self.tot_num_x_type[:i]) for i in xrange(len(self.tot_num_x_type))] # indexes of the neurons divided per type
		self.halves = [w+(u-w)/2 for (u,w) in zip(self.num_types[1:],self.num_types[0:-1])] # half number of cells
		self.vect_index = [range(i,j) for (i,j) in zip(self.num_types[0:-1],self.num_types[1:])] # vector of indexes 
		self.left_index = [range(i,j) for (i,j) in zip(self.num_types[0:-1],self.halves)] # left indexes
		self.right_index = [range(i,j) for (i,j) in zip(self.halves,self.num_types[1:])] # right indexes
		self.colors = ["#ffd232","#ff0000","#4646b4","#00aadc","#96501e","#00963c","#ffaa8c"] # colors

		# ===== NUMERICAL INTEGRATION PARAMETERS =====
		self.varDt = False # numbrical integration: True=variable time step, False=fixed time step (standard is 0.01)
		self.atol = 1e-5 # absolute tolerance if varDt=True
		self.rtol = 1e-5 # relative tolerance if varDt=True
		self.dt = 0.01 # time step if varDt=False
		self.RecAll = 0 # 0="no variable recorded" - 1="record only time-voltage" - 2:"record also synaptic currents"
		self.time_end = 1000 # total simulation time 

		# ===== PARAMETERS FOR INJECTING CURRENT TO SELECTED NEURONS =====
		self.pos_active = 30 # index of RBs that are activated using a pulse of current to start the swimmming dynamics
		self.n_active = 2 # number of active RBs - i.e. number of RBs receiving the current pulse
		self.delay = 50.0 # ms from the start of current injection
		self.duration = 5.0 # ms of inj duration
		self.amplitude_mean = 0.15 # in nA
		self.amplitude_std = 0.0 # in nA

		# ===== SYNAPTIC-ANATOMICAL PARAMETERS =====
		self.sa_prop = 0.02 # variability in the cell RC position 
		self.fixed_delay = 1.0 # fixed synaptic delay 
		self.var_delay = 0.0035 # variability in the synaptic delay to multiply by the RC pre-post synaptic cell RC distance
		self.std_on = 1

		# ===== SYNAPTIC PARAMETERS =====
		ampa_strength = 0.593e-3 # standard AMPA strengths 
		inh_strength = 0.435e-3 # standard inh strengths 
		nmda_strength = 0.29e-3 # standard NMDA strengths

		'''
		w = synaptic parameters given in the following format "pre neuron -> post" in each n-pla (...): 
		[syn_type, mean strength, std strength, delay ("distance"=RC distant dependentif, otherwise you can set it to a fixed delay)]
		'''
		self.w = { 
			# pre rb
			"rb -> dla": (["ampa", 8.0e-3, 0.05*8.0e-3, "distance"],None),
			"rb -> dlc": (["ampa", 8.0e-3, 0.05*8.0e-3, "distance"],None),
			"rb -> din": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),

			# pre din
			"din -> din": (["ampa", ampa_strength, ampa_strength*0.05, "distance"], ["nmda", 0.15e-3, 0.0075*1e-3, "distance"]),
			"din -> cin": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"din -> dlc": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"din -> mn": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"din -> dla": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"din -> ain": (["ampa", 0.1e-3, 0.05*0.1*1e-3, "distance"],None),

			# pre ain
			"ain -> dlc": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"ain -> ain": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"ain -> cin": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"ain -> din": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"ain -> mn": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"ain -> dla": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
	
			# pre cin
			"cin -> dlc": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"cin -> ain": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"cin -> cin": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"cin -> din": (["inh", inh_strength, inh_strength*0.05, "distance"],None), 
			"cin -> mn": (["inh", inh_strength, inh_strength*0.05, "distance"],None),
			"cin -> dla": (["inh", inh_strength, inh_strength*0.05, "distance"],None),

			# pre mn
			"mn -> ain": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"mn -> cin": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"mn -> din": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),
			"mn -> mn": (["ampa", ampa_strength, ampa_strength*0.05, "distance"],None),

			# start as in Roberts et al, 2014
			"dlc -> din": (None,["nmda", 1.0e-3, 0.05*1.0e-3, "distance"]),
			"dlc -> ain": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dlc -> cin": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dlc -> mn": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dlc -> dlc": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),

			"dla -> dla": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dla -> din": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dla -> ain": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dla -> cin": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dla -> mn": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			"dla -> dlc": (None,["nmda", nmda_strength, nmda_strength*0.05, "distance"]),
			}
	

def create_params():
	params = SimulationParams()
	return params






