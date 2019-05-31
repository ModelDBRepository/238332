
'''
This file contains utility rountines are used in the main.py file. These routines include
plotting tools for visualization, generation of connectivity and current clamp.
'''

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rnd
import random

# import parameters from param file
import param
par = param.create_params()
path = par.path
cell_types = par.cell_types
num_types = par.num_types
halves = par.halves
vect_index = par.vect_index
left_index = par.left_index
right_index = par.right_index
colors = par.colors

# return spike train of a simulated cell 
def SpikeTrain(cellist):
	return [np.array(cell.record['spk']) for cell in cellist]

# return time and voltage of cell
def TimeVoltageTrace(cell):
	return (np.array(cell.record['t']),np.array(cell.record['vm']))

# plot single spike train
def plot_spk_train(pos,st,I='',colors=''):
	for (i,spike_train) in enumerate(st):
		plt.plot(spike_train,pos[i]*np.ones_like(spike_train),marker='.',markersize=5,markerfacecolor=colors,markeredgecolor=colors,linestyle='None')
	if len(I):
		plt.xlim(I)

# current clamp isntance
def inj_current(cells,delay,dur,amp_mean=0,amp_std=0):
	for cell in cells:
		cell.CC.delay = delay
		cell.CC.dur = dur
		cell.CC.amp = rnd.normal(amp_mean,amp_std)


# spike train raster plot 
def plotLeftRightSpikeTrain(cellist,I):
	st=np.array([train.tolist() for train in SpikeTrain(cellist)])
	plt.subplot(2,1,2)
	for i in xrange(len(left_index)):
		if i!=15: # do not plot ith cell type
			plot_spk_train([cell.pos for cell in [cellist[x] for x in left_index[i]]],st[left_index[i]],I,colors=colors[i])
			plt.ylim([500,2000])
			plt.ylabel("Left")
	plt.subplot(2,1,1)
	for i in xrange(len(right_index)):
		if i!=15: # do not plot ith cell type
			plot_spk_train([cell.pos for cell in [cellist[x] for x in right_index[i]]],st[right_index[i]],I,colors=colors[i])
			plt.ylim([500,2000])
			plt.ylabel("Right")

# plot voltage traces of selected cellist vector with offset voltage
def plotLeftRightVoltageOffset(cellist,I,offset=0.0):
	n = len(cellist)
	plt.figure()
	for i in xrange(n):
		if i<n/2:
			plt.subplot(2,1,2)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v+y_offset*i,color=cellist[i].color,linewidth=1.0)
			plt.xlim(I)
			plt.ylabel("Left")
		else:
			plt.subplot(2,1,1)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v+y_offset*i,color=cellist[i].color,linewidth=1.0)
			plt.xlim(I)
			plt.ylabel("Right")
	plt.subplot(2,1,2)
	plt.ylim([-80,40+offset*n])
	plt.subplot(2,1,2)
	plt.ylim([-80,40+offset*n])

# plotting voltage traces of selected cellist vector
def plotLeftRightVoltage(cellist,idx,I):
	n = len(cellist)
	#plt.figure()
	for i in idx:
		if cellist[i].body_side==1:
			plt.subplot(4,1,4)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
			plt.xlim(I)
			plt.ylabel("Left")
		elif cellist[i].body_side==2:
			plt.subplot(4,1,1)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
			plt.xlim(I)
			plt.ylabel("Right")
		else:
			print "Error no body side declaration"

	st=np.array([train.tolist() for train in SpikeTrain(cellist)])
	plt.subplot(4,1,3)
	for i in xrange(len(left_index)):
		plot_spk_train([cell.pos for cell in [cellist[x] for x in left_index[i]]],st[left_index[i]],I,colors=colors[i])
	plt.subplot(4,1,3)
	plt.ylim([500,2000])
	plt.ylabel("Left")

	plt.subplot(4,1,2)
	for i in xrange(len(right_index)):
		plot_spk_train([cell.pos for cell in [cellist[x] for x in right_index[i]]],st[right_index[i]],I,colors=colors[i])
	plt.subplot(4,1,2)
	plt.ylim([500,2000])
	plt.ylabel("Right")

	plt.subplot(4,1,4)
	plt.ylim([-80,40])
	plt.subplot(4,1,1)
	plt.ylim([-80,40])


# destroy all cellists - calling method contained in the cellist class
def destroy(cellist):
	for cell in cellist:
		cell.destroy()











