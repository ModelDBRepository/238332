
import numpy as np
import scipy.stats as stats
import scipy.stats.mstats as mstats
import matplotlib.pylab as plt
plt.ion()

path="/home/andrea/Desktop/old codes connectome/neuron tadpole original connectome old/figures"
n=200


num_types=[0,126,230,366,750,986,1324,1382] # number of cells per cell type
halves=[w+(u-w)/2 for (u,w) in zip(num_types[1:],num_types[0:-1])]
cell_types=["rb","dlc","ain","cin","din","mn","dla"]
vect_index=[range(i,j) for (i,j) in zip(num_types[0:-1],num_types[1:])]
left_index=[range(i,j) for (i,j) in zip(num_types[0:-1],halves)]
right_index=[range(i,j) for (i,j) in zip(halves,num_types[1:])]
colors=["#ffd232","#ff0000","#4646b4","#00aadc","#96501e","#00963c","#ffaa8c"]

bound1=45
bound2=75

noise=" new"
start=1

pos=np.loadtxt(path+"/../data/pos.txt")
pos=np.array(pos[left_index[3]])

mn_mean_t_ana=[]
indeg_cins=[]
for num in range(start,n+1):
	print num
	spk=np.load(path+"/anatomical 1.5 sec"+noise+"/spk_ana"+str(num)+".npy")
	mn_t_ana=[]
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_ana.append(tmp)
	mn_mean_t_ana.append(np.median(mn_t_ana))

	A=np.load(path+"/anatomical 1.5 sec"+noise+"/A_anatomical"+str(num)+".npy")
	indeg_cins.append(sum(A[np.ix_(left_index[4],left_index[3])]))

plt.figure(1)
mean_indeg=np.mean(indeg_cins,axis=0)
std_indeg=np.std(indeg_cins,axis=0)
plt.subplot(1,2,1)
plt.fill_between(pos,mean_indeg-std_indeg,mean_indeg+std_indeg,alpha=0.5)
plt.plot(pos,mean_indeg,'k',lw=3)
plt.ylabel("in-degree from dINs")
plt.xlim([500,2000])
plt.ylim([3,30])
plt.axhline(13, color='r', lw=2, alpha=1)

mn_mean_t_prob=[]
indeg_cins=[]
for num in range(start,n+1):
	print num
	spk=np.load(path+"/probabilistic 1.5 sec"+noise+"/spk_prob"+str(num)+".npy")
	mn_t_prob=[]	
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_prob.append(tmp)
	mn_mean_t_prob.append(np.median(mn_t_prob))

	A=np.load(path+"/probabilistic 1.5 sec"+noise+"/A_probabilistic"+str(num)+".npy")
	indeg_cins.append(sum(A[np.ix_(left_index[4],left_index[3])]))

plt.figure(1)
plt.subplot(1,2,2)
plt.fill_between(pos,mean_indeg-std_indeg,mean_indeg+std_indeg,alpha=0.5,label="mean+std")
plt.plot(pos,mean_indeg,'k',lw=3)
plt.xlim([500,2000])
plt.ylim([3,30])
plt.axhline(13, color='r', lw=2, alpha=1)



plt.figure(2)
col=["grey","k"]
mews=2
size=15
frame = plt.gca()
dx=40
n,bins,patches=plt.hist([mn_mean_t_ana,mn_mean_t_prob],dx,label=['Anatomical','Probabilistic'],color=[col[0],col[1]])
plt.ylabel("Number of Simulations")
plt.xlabel("Swimming Period")
#frame.axes.get_xaxis().set_visible(False)
plt.xlim([50,72])
#plt.ylim([0,np.max(n)+1])
plt.legend(loc=2)









