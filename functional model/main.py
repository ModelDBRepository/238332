
'''
This file is the main file that simulates the functional model of spiking neurons of the swimming spinal 
network in tadpoles. There are several routines and subroutines, that briefly summarize in creation of cells and 
their connectivity, launching simulations, plotting results, parallelize simulations across a 
selected number of cores. The 
'''

from nrnTemplate.CellTypes.CellTemplateCommon import dINcell, MNcell
import matplotlib.pylab as plt
import numpy as np
import random
import numpy.random as rnd
import mp_util
from util import *
from neuron import h

# create folder named by the current time where to save simulation results 
import time, datetime, os
from shutil import copyfile
path_tmp = "figures/"
today = datetime.datetime.now()
todaystr = today.isoformat()
os.mkdir(path_tmp+todaystr)
save_path = path_tmp+todaystr+"/"

# initialize parameter from param.py 
import param
par = param.create_params()
Q = par.P
posi = par.pos
path = par.path
RecAll = par.RecAll
varDt = par.varDt
atol = par.atol
rtol = par.rtol
dt = par.dt
time_end = par.time_end
cell_types = par.cell_types
num_types = par.num_types
halves = par.halves
vect_index = par.vect_index
left_index = par.left_index
right_index = par.right_index
colors = par.colors
fixed_delay = par.fixed_delay
sa_prop = par.sa_prop
var_delay = par.var_delay
std_on = par.std_on
w = par.w
n_active = par.n_active
pos_active = par.pos_active
delay = par.delay
duration = par.duration
amplitude_mean = par.amplitude_mean
amplitude_std = par.amplitude_std


# Cells Creation 
def CellCreation(RecAll=0,varDt=False,sim_num=1,sim_type=1):
	if sim_type==0:
		pos = posi	
	else:
		pos = np.load(path+"pos"+str(sim_num)+".npy").tolist()
	cellist=[]
	for i in xrange(len(cell_types)):
		for j in vect_index[i]:
			if cell_types[i] == "din":
				cellist.append(dINcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol,theta=0.0))
			elif cell_types[i] == "cin":
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			else:
				cellist.append(MNcell(RecAll=RecAll,varDt=varDt,atol=atol,rtol=rtol))

			cellist[-1].whatami=cell_types[i]
			cellist[-1].color=colors[i]
			cellist[-1].index=len(cellist)-1
			cellist[-1].type_num=i
			cellist[-1].pos=rnd.normal(pos[j],abs(pos[j])*sa_prop)
			if j in left_index[i]:
				cellist[-1].body_side=1
			else:
				cellist[-1].body_side=2
	return cellist


# create synapse between neural types
def connection(pre,post,dist):
	key = pre.whatami + " -> " + post.whatami
	if w.has_key(key):
		specs = w[key]
		for spec in specs:
			if spec != None:
				syn_type = spec[0]
				w_mean = spec[1]
				w_std = spec[2]
				if w_std != 0.0:
					weight = rnd.normal(w_mean,w_std*std_on)
				else:
					weight = w_mean

				if spec[3]=="distance":
					distance = fixed_delay+dist*var_delay
				else:
					distance = spec[3]
				pre.connect(post, syn_type, w=weight, delay=distance)
				
	else:
		print "Something wrong with this connection type\n"
		print key
		raise Exception("Connection created is not included in the list of possible connections")

# create gap junctions between dINs
def GapJunction(cellist):
	for i in left_index[4]:  
		for j in left_index[4]:
			if i!=j and abs(cellist[i].pos-cellist[j].pos)<100.0: 
				cellist[i].connect(cellist[j],"gap",gmax=0.2e-3)

	for i in right_index[4]:
		for j in right_index[4]:
			if i!=j and abs(cellist[i].pos-cellist[j].pos)<100.0: 
				cellist[i].connect(cellist[j],"gap",gmax=0.2e-3)

# no Commissural cINs
def ConfigNoCommissural(cellist,P):
	n=len(P)
	for i in xrange(n):
		for j in xrange(n):
			if cellist[i].whatami=="cin":
				P[i,j]=0
	return P

# cIN-dIN subnetwork - remove all mn and aIN connections 
def ConfigSubnet_dIN_cIN(cellist,P):
	n=len(P)
	for i in xrange(n):
		for j in xrange(n):
			if cellist[j].whatami=="mn" or cellist[j].whatami=="ain":
				P[i,j]=0
	return P

# no ascending dINs
def ConfigNoAscendingdINs(cellist,P):
	n=len(P)
	for i in xrange(n):
		for j in xrange(n):
			if cellist[i].whatami=="din" and cellist[j].whatami=="din" and cellist[i].pos>cellist[j].pos:
				P[i,j]=0
	return P

# probabilistic connectome
def CreateConfigAdj(cellist,sim_num,sim_type=1):
	if sim_type==0: # anatomical connectome
		A = np.load(path+"A_connectome"+str(sim_num)+".npy")
	else:
		if sim_type==1: # probabilistic connectome
			P = Q
		elif sim_type==2: # 2=cIN-dIN probabilistic subnetwork 
			P = Q
			P = ConfigSubnet_dIN_cIN(cellist,P)
		elif sim_type==3: # 3=no commissural connections
			P = Q
			P = ConfigNoCommissural(cellist,P)
		elif sim_type==4: # 4=no ascending dINs
			P = Q
			P = ConfigNoAscendingdINs(cellist,P)

		n=len(P)
		A = np.zeros((n,n))
		for i in xrange(n):
			for j in xrange(n):
				if P[i,j]>rnd.rand():
					A[i,j]=1

	n=len(cellist)
	for i in xrange(n):
		for j in xrange(n):
			if A[i,j]:
				connection(cellist[i],cellist[j],abs(cellist[i].pos-cellist[j].pos))
	np.save(save_path+"A"+str(sim_num)+".npy",A)

# run one simulation of the swimming network - RB neurons are initialized before the simulation starts
def SwimmingSimulation(tstop,sim_num,sim_type=1):
	t_start = time.time()

	seed = random.getrandbits(32) # set a fixed seed if desired
	rnd.seed(seed)
	with open(save_path+"seed", "wt") as f_seed:
		f_seed.write(str(seed))
	print "Set Random Seed = {0}.".format(seed)

	print "Running Simulation "+str(sim_num)+"..."
	
	print "Create Cells ..."
	cellist=CellCreation(RecAll=RecAll,varDt=varDt,sim_num=sim_num,sim_type=sim_type)
	print "Cells Created"

	print "Create Connectivity ..."
	CreateConfigAdj(cellist,sim_num)
	GapJunction(cellist)
	print "Connectivity Created"

	print "Run Numerical Integration..."
	inj_current(cellist[pos_active:pos_active+n_active],delay,duration,amp_mean=amplitude_mean,amp_std=amplitude_std) 
	RunSim(tstop=tstop,dt=dt)
	print "End of the Integration"

	# plot voltage of selected neurons and spike train divided between left-right sides
	plt.figure()
	plotLeftRightVoltage(cellist,[left_index[4][59],left_index[3][86],right_index[4][59],right_index[3][86]],[0,tstop])

	# save spike trains
	np.save(save_path+"spikes"+str(sim_num)+".npy",[np.array(cell.record["spk"]) for cell in cellist])

	t_end = time.time()
	print "Simulation Took {0}s.".format(t_end-t_start)

# create simulation routine that depends on the simulation number (used by RunManySim() routine)
def SimForAll(sim_num):
	output=SwimmingSimulation(time_end,sim_num,sim_type=1)
	return output

# run multiple parallel simulations using mp_util.py with selected number of cores (num_process)
def RunManySim(num_process):
	# save files param.py, util.py and main.py to record simulation parameters and routines
	file_main = "param.py"
	if ~os.path.isfile(save_path+file_main):
		copyfile(file_main,save_path+file_main)
	file_util = "util.py"
	if ~os.path.isfile(save_path+file_util):
		copyfile(file_util,save_path+file_util)
	file_param = "main.py"
	if ~os.path.isfile(save_path+file_param):
		copyfile(file_param,save_path+file_param)
	processes = mp_util.UniqueProcessMap(num_process)
	I = range(1,201)
	out=processes.map(SimForAll,I)
	return out

# run simulation using Euler/CVode method
def RunSim(v_init=-80.0,tstop=0.0,dt=0.01):
	t_start = time.time()
	h.dt = dt
	h.t = 0.0
	h.finitialize(v_init)
	while h.t<tstop:
		h.fadvance()








