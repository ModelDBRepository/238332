
from abc import ABCMeta
import numpy as np
from neuron import h

class Cell(object):
	__metaclass__ = ABCMeta
	
	def __init__(self, cap=4.0e-6, RecAll=0, varDt=False, atol=1e-5, rtol=1e-5, theta=0.0, channels={}):
		'''
			This class defines a single compartment cell standard class with the following parameters
			- cap: capacitance
			- RecAll: 0 save nothing, 1 saves time and voltage, 2 saves also synaptic currents
			- varDt: True set the use of local time step integrator, False uses predefined step dt
			- atol, rtol: in case varDt=1 sets the tolerances of the CVode numerical integration
			- theta: threshold for spike transmission
			- channels: active channels
		'''
	
		# position and side
		self.pos = 0.0
		self.side = 0 # 0=unassigned, 1=left, 2=right
		self.index = None

		# param
		self.RecAll = RecAll
		self.theta = theta
	
		# soma anatomy
		self.soma=h.Section(cell=self)
		self.soma.nseg=1
		self.soma.L=100.0 # um (?)
		self.surface_area=1.0 # cm^2
		self.soma.diam=(self.surface_area*1e8)/(np.pi*self.soma.L)
		self.soma.cm=cap

		# soma ionic channels
		for ch in channels:
			self.soma.insert(ch)
		
		# set current clamp and voltage clamp 
		self.CC = h.IClamp(self.soma(0.5))
		self.VC = h.VClamp(self.soma(0.5))

		# record spike times
		self.record = {}
		self.threshold = 0.0
		self.nc_spike = h.NetCon(self.soma(0.5)._ref_v,None,self.threshold,0.0,0.0,sec=self.soma)
		self.record['spk'] = h.Vector()
		self.nc_spike.record(self.record['spk'])

		# vector of connections
		self.connlist = []
		self.gap_list=[]

		if varDt:
			# variable time step integrator
			self.Hines = h.CVode()
			self.Hines.active(1)
			self.Hines.use_local_dt(1)
			self.Hines.atol(atol)
			self.Hines.rtol(rtol)

		if RecAll>0:
			self.record['t']=h.Vector()
			self.record['vm']=h.Vector()
			if varDt:
				self.Hines.record(self.soma(0.5)._ref_v,self.record['vm'],self.record['t'],sec=self.soma)
				if RecAll==2:
					#self.record['i_ca']=h.Vector()
					self.record['i_kf']=h.Vector()
					#self.record['i_ca'].record(self.soma(0.5)._ref_ica)
					self.record['i_kf'].record(self.soma(0.5)._ref_ik)
			else: 
				self.record['t'].record(h._ref_t)
				self.record['vm'].record(self.soma(0.5)._ref_v)
				if RecAll==2:
					#self.record['i_ca']=h.Vector()
					self.record['i_kf']=h.Vector()
					#self.record['i_ca'].record(self.soma(0.5)._ref_ca)
					self.record['i_kf'].record(self.soma(0.5)._ref_ik)


	def connect(self,dest,InputSynName,w=0.0,delay=0.0,gmax=0.0):
		if InputSynName.rfind('ampa')>=0:
			if hasattr(dest,"syn_ampa")==0:
				dest.syn_ampa = h.syn_ampa(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa,self.theta,delay,w,sec=self.soma))
			if self.RecAll==2:
				dest.record['i_ampa']=h.Vector() 
				dest.record['i_ampa'].record(dest.syn_ampa._ref_i,sec=self.soma)
		if InputSynName.rfind('ampa_dale')>=0:
			if hasattr(dest,"syn_ampa_dale")==0:
				dest.syn_ampa_dale = h.syn_ampa_dale(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa_dale,self.theta,delay,w,sec=self.soma))
			if self.RecAll==2:
				dest.record['i_ampa_dale']=h.Vector() 
				dest.record['i_ampa_dale'].record(dest.syn_ampa_dale._ref_i,sec=self.soma)
		if InputSynName.rfind('nmda')>=0:
			if hasattr(dest,"syn_nmda")==0:
				dest.syn_nmda = h.syn_nmda(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_nmda,self.theta,delay,w,sec=self.soma))
			if self.RecAll==2:
				dest.record['i_nmda']=h.Vector() 
				dest.record['i_nmda'].record(dest.syn_nmda._ref_i,sec=self.soma)
		if InputSynName.rfind('inh')>=0:
			if hasattr(dest,"syn_inh")==0:
				dest.syn_inh = h.syn_inh(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh,self.theta,delay,w,sec=self.soma))
			if self.RecAll==2:
				dest.record['i_inh']=h.Vector()
				dest.record['i_inh'].record(dest.syn_inh._ref_i,sec=self.soma)
		if InputSynName.rfind('inh_mhr')>=0:
			if hasattr(dest,"syn_inh_mhr")==0:
				dest.syn_inh_mhr = h.syn_inh_mhr(dest.soma(0.5))
				dest.syn_inh_mhr.tau_c = 15.0
			dest.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh_mhr,self.theta,delay,w,sec=self.soma))
			if self.RecAll==2:
				dest.record['i_inh']=h.Vector()
				dest.record['i_inh'].record(dest.syn_inh._ref_i,sec=self.soma)
		elif InputSynName.rfind('gap')>=0:
			dest.gap_list.append(h.Gap(dest.soma(0.5)))
			dest.gap_list[-1].gmax=gmax
			dest.connlist.append(h.setpointer(self.soma(0.5)._ref_v,'vgap',dest.gap_list[-1]))
			if self.RecAll==2:
				dest.record['i_gap']=h.Vector() 
				dest.record['i_gap'].record(dest.gap._ref_i,sec=self.soma)

	# destroy
	def destroy(self):
		del self.record
		del self.soma


class dINcell(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"
		
		# channel properties and cell init
		self.channels = {"pas","dIN_na","dIN_kFast","dIN_kSlow","dIN_ca"}

		# init cell		
		Cell.__init__(self, cap=10.0e-6, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -52.0, 
		                   "gmax_lk": 1.405e-9,
		                   "erev_na": 50.0,
		                   "gmax_na": 240.5e-9, # Down to 210.5e-9 reduces mid-cycle dINs
		                   "erev_k": -81.5,
		                   "gmax_kf": 12.0e-9,
		                   "gmax_ks": 9.6e-9,
		                   "T_ca": 300.0,
		                   "Sin_ca": 100e-9,
		                   "Sout_ca": 10e-6,
		                   "perm_ca": 1.425e-10} # 0.725e-10 gives a decent reduction in mid-cycle dINs

		# ion channels
		self.soma.L=100
		self.soma.cm=10.0e-6
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_dIN_na=self.parameters["gmax_na"]
		self.soma.gmax_dIN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_dIN_kSlow=self.parameters["gmax_ks"]
		self.soma.T_dIN_ca=self.parameters["T_ca"]
		self.soma.perm_dIN_ca=self.parameters["perm_ca"]
		self.soma.Sin_dIN_ca=self.parameters["Sin_ca"]
		self.soma_Sout_dIN_ca=self.parameters["Sout_ca"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]



class MNcell(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="mn"
		self.color="#00963c"
		
		# channel properties and cell init
		self.channels = {"pas","MN_na","MN_kFast","MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=10.0e-6, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"L": 100.0,
		                  "erev_lk": -61.0, 
		                  "gmax_lk": 2.4691e-9, 
		                  "erev_na": 50.0,
		                  "gmax_na": 110e-9,
		                  "erev_k": -80.0,
		                  "gmax_kf": 8e-9,
		                  "gmax_ks": 1e-9}

		# ion channels
		self.soma.L=100
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_MN_na=self.parameters["gmax_na"]
		self.soma.gmax_MN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_MN_kSlow=self.parameters["gmax_ks"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]











		
