
import numpy as np
import scipy.stats as stats
import scipy.stats.mstats as mstats
import matplotlib.pylab as plt
from scipy import stats
plt.ion()

# ===== LINEAR REGRESSION MODEL =====
n=100
path_anat="/home/andrea/Desktop/old codes connectome/neuron tadpole original connectome old/figures/anatomical 1.5 sec new" # path of the anatomical connectomes 
path_prob="/home/andrea/Desktop/old codes connectome/neuron tadpole original connectome old/figures/probabilistic 1.5 sec new" # path of the probabilistic connectomes

num_types=[0,126,230,366,750,986,1324,1382] # number of cells per cell type
halves=[w+(u-w)/2 for (u,w) in zip(num_types[1:],num_types[0:-1])]
cell_types=["rb","dlc","ain","cin","din","mn","dla"]
vect_index=[range(i,j) for (i,j) in zip(num_types[0:-1],num_types[1:])]
left_index=[range(i,j) for (i,j) in zip(num_types[0:-1],halves)]
right_index=[range(i,j) for (i,j) in zip(halves,num_types[1:])]
colors=["#ffd232","#ff0000","#4646b4","#00aadc","#96501e","#00963c","#ffaa8c"]

bound1=50
bound2=75
thresh1=0.0
thresh2=13.0
k=2
r=1
noise=" new"
start=1

mn_mean_t_ana=[]
din_mean_ana=[]
cin_mean_ana=[]

for num in range(start,n+1):
	spk=np.load(path_anat+"/spk_ana"+str(num)+".npy")
	mn_t_ana=[]	
	count_mn=0
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_ana.append(tmp)
	mn_mean_t_ana.append(np.median(mn_t_ana))

	A=np.load(path_anat+"/A_anatomical"+str(num)+".npy")
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[4])])>=thresh1)
	conn=sum(A[np.ix_(750+idx[0],vect_index[4])])
	din_mean_ana.append(np.mean(conn))
	
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[3])])>=thresh2)
	conn=sum(A[np.ix_(366+idx[0],vect_index[4])])
	cin_mean_ana.append(np.mean(conn))

mn_mean_t_ana=np.array(mn_mean_t_ana)
din_mean_ana=np.array(din_mean_ana)
cin_mean_ana=np.array(cin_mean_ana)


din_mean_prob=[]
cin_mean_prob=[]
mn_mean_t_prob=[]

for num in range(start,n+1):
	spk=np.load(path_prob+"/spk_prob"+str(num)+".npy")
	mn_t_prob=[]	
	count_mn=0
	for i in vect_index[5]:
		if len(spk[i])>20	:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_prob.append(tmp)
	mn_mean_t_prob.append(np.median(mn_t_prob))	

	A=np.load(path_prob+"/A_probabilistic"+str(num)+".npy")
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[4])])>=thresh1)
	conn=sum(A[np.ix_(750+idx[0],vect_index[4])])
	din_mean_prob.append(np.mean(conn))

	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[3])])>=thresh2)
	conn=sum(A[np.ix_(366+idx[0],vect_index[4])])
	cin_mean_prob.append(np.mean(conn))

	
mn_mean_t_prob=np.array(mn_mean_t_prob)
din_mean_prob=np.array(din_mean_prob)
cin_mean_prob=np.array(cin_mean_prob)

col=["grey","k"]
mews=2
size=15


s=np.concatenate([mn_mean_t_ana,mn_mean_t_prob])
r1=np.concatenate([din_mean_ana,din_mean_prob])
r3=np.concatenate([cin_mean_ana,cin_mean_prob])

x=[r1,r3]
y=s

X = np.column_stack(x+[[1]*len(x[0])])
beta_hat = np.linalg.lstsq(X,y)[0]
xnew=np.dot(X,beta_hat)
slope, intercept, r_value, p_value, std_err = stats.linregress(xnew,y)

# ===== PREDICTION TEST =====
noise=" new"
start=100
n=200
true=[]
predict=[]
r=beta_hat

for num in range(start,n+1):
	spk=np.load(path_anat+"/spk_ana"+str(num)+".npy")
	mn_t_ana=[]	
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_ana.append(tmp)
	true.append(np.median(mn_t_ana))

	A=np.load(path_anat+"/A_anatomical"+str(num)+".npy")
	
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[4])])>=thresh1)
	conn=sum(A[np.ix_(750+idx[0],vect_index[4])])
	din_mean_ana=np.mean(conn)
	
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[3])])>=thresh2)
	conn=sum(A[np.ix_(366+idx[0],vect_index[4])])
	cin_mean_ana=np.mean(conn)
	
	x1=din_mean_ana
	x2=cin_mean_ana
	predict.append(x1*r[0]+x2*r[1]+r[2])


for num in range(start,n+1):
	spk=np.load(path_prob+"/spk_prob"+str(num)+".npy")
	mn_t_prob=[]	
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				mn_t_prob.append(tmp)
	true.append(np.median(mn_t_prob))

	A=np.load(path_prob+"/A_probabilistic"+str(num)+".npy")
	
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[4])])>=thresh1)
	conn=sum(A[np.ix_(750+idx[0],vect_index[4])])
	din_mean_prob=np.mean(conn)
	
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[3])])>=thresh2)
	conn=sum(A[np.ix_(366+idx[0],vect_index[4])])
	cin_mean_prob=np.mean(conn)
	
	x1=din_mean_prob
	x2=cin_mean_prob

	predict.append(x1*r[0]+x2*r[1]+r[2])


plt.figure()
x=np.array(predict)
y=np.array(true)


idx=n-start+1
plt.plot(x[xrange(0,idx)],y[xrange(0,idx)],'.',color=col[0],ms=15,markeredgecolor=col[0],mew=mews,label="Anatomical")
plt.plot(x[xrange(idx,2*idx)],y[xrange(idx,2*idx)],'.',color=col[1],ms=15,markeredgecolor=col[1],mew=mews,label="Probabilistic")
plt.plot([min(x),max(x)],slope*np.array([min(x),max(x)])+intercept,color="b",ms=size,markeredgecolor=col[1],mew=mews,label="Perfect Match")
plt.xlabel("Predicted Period")
plt.ylabel("True Period")
plt.ylim([53,71])
plt.xlim([53,71])

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print stats.pearsonr(x,y)


















