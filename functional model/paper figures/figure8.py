
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import scipy.stats as stats
import matplotlib.transforms as mtransforms

num_types=[0,126,230,366,750,986,1324,1382] # number of cells per cell type
halves=[w+(u-w)/2 for (u,w) in zip(num_types[1:],num_types[0:-1])]
cell_types=["rb","dlc","ain","cin","din","mn","dla"]
vect_index=[range(i,j) for (i,j) in zip(num_types[0:-1],num_types[1:])]
left_index=[range(i,j) for (i,j) in zip(num_types[0:-1],halves)]
right_index=[range(i,j) for (i,j) in zip(halves,num_types[1:])]

pos=np.loadtxt('pos_1000.txt')
P=np.loadtxt("P_1000.txt")
P=P.transpose()

def firing_rate(st,interval):
	return [len([item for item in spikes if item<interval[1] and item>interval[0]]) for spikes in st]

idx_cin=np.array(vect_index[3])
tstop=2000
path="/home/andrea/Desktop/old codes connectome/neuron tadpole original connectome old/figures/probabilistic 1.5 sec new/"
#path='/home/andrea/Desktop/computer_right/PyNN_Tadpole_Simulation (copy)/figures prob model/norm prob 6 sec new'
I=[1000,2000]

pos_bin=np.linspace(500,2000,50)
cin_pos=pos[vect_index[3]]
cin_pos=np.array(cin_pos)

unr_cin_all=[]
for i in range(1,101):
	st=np.load(path+"spk_prob"+str(i)+".npy")
	cin=st[vect_index[3]]
	cin_fr=np.array(firing_rate(cin,I))
	tmp=[]
	for k in cin_fr:
		if k<10:
			tmp.append(1)
		else:
			tmp.append(0)
	unr_cin_all.append(tmp)


# ===== PLOTS =====
fig, (ax1, ax2) = plt.subplots(2,1)
inc_conn_din=P[np.ix_(vect_index[4],vect_index[3])].sum(axis=0)
n=len(inc_conn_din)
grp1=[i for i in xrange(n) if inc_conn_din[i]<=15]
grp2=[i for i in xrange(n) if inc_conn_din[i]>15]

ax1.plot(cin_pos[grp1],inc_conn_din[grp1],marker='+',linestyle="None",color="r",mew=5, ms=1)
ax1.plot(cin_pos[grp2],inc_conn_din[grp2],marker='+',linestyle="None",color="#00aadc",mew=5, ms=1)

ax1.axhline(13, color='k', lw=1, alpha=0.1)
ax1.axhline(15, color='k', lw=1, alpha=0.1)
trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)

x=np.linspace(500,2000,1e3)
ax1.fill_between(x, 13, 15, facecolor='green', alpha=0.2)
ax1.set_xlim([500,2000])


plt.figure(1)
plt.subplot(2,1,2)
plt.plot(pos[idx_cin],1-np.mean(unr_cin_all,axis=0),marker=".",markersize=5,linestyle='None',mew=3, ms=3)
plt.xlim([500,2000])
plt.ylim([0,1])

size=15
mews=2
plt.figure(2)
xx=inc_conn_din
yy=1-np.mean(unr_cin_all,axis=0)
plt.plot(xx,yy,marker=".",markersize=5, markeredgecolor="k", linestyle='None',mew=5, ms=5)
slope, intercept, r_value, p_value, std_err = stats.linregress(xx,yy)
plt.plot([min(xx),max(xx)],slope*np.array([min(xx),max(xx)])+intercept,color="b",ms=size,lw=2,markeredgecolor="r",mew=mews,label="Perfect Match")
plt.xlabel("# of dIN-cIN connections")
plt.ylabel("Unreliable cINs")
plt.ylim([min(yy)-0.01,max(yy)+0.01])
plt.xlim([min(xx)-0.2,max(xx)+0.2])
















