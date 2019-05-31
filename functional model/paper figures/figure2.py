
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

path="../probabilistic model/anatomical adjacency matrixes/"
num=1382
num_x_type=[0,126, 104, 136, 384, 236, 338, 58]
halves=np.divide(num_x_type,2)
types=['RB','DLC','aIN','cIN','dIN','MN','DLA']
order_best=[0,6,1,2,3,4,5]

grid_color1="#87CEFA"
grid_color2="#87CEFA"
lw=0.5
tr=0.5

P=np.load(path+"P.npy")
[n,m]=P.shape
Q=np.zeros(P.shape)

v_ind=[]
for i in xrange(len(num_x_type)+1):
	if i!=0:
		v_ind.append(np.sum(num_x_type[0:i]))

new_idx=[]
for i in xrange(len(order_best)):
	new_idx.append(range(v_ind[order_best[i]],v_ind[order_best[i]+1]))

new_idx=np.concatenate((new_idx)).tolist()
Q=P[new_idx,:]
Q=Q[:,new_idx]

nch=1
colors=["#ffd232","#ffaa8c","#ff0000","#4646b4","#00aadc","#96501e","#00963c"]
num_x_type=[0, 126, 58, 104, 136, 384, 236, 338]
halves=np.divide(num_x_type,2)
types=['RB','dla','dlc','aIN','cIN','dIN','mn']

plt.matshow(Q.transpose(),cmap='Greys_r')

plt.plot([-100,num-0.5],[num-0.5,num-0.5],color=grid_color1)
plt.plot([-100,num-0.5],[-100,-100],color=grid_color1)
plt.plot([num-0.5,num-0.5],[-100,num-0.5],color=grid_color1)
plt.plot([-100,-100],[-100,num-0.5],color=grid_color1)
plt.text(-60,-40,'P', fontsize=20)

for i in xrange(len(num_x_type)):
	line_half=np.sum(num_x_type[0:i])+halves[i]-tr
	line=np.sum(num_x_type[0:i])-tr
	plt.plot([0,num-tr],[line_half,line_half],':',color=grid_color2,linewidth=lw)
	plt.plot([line_half,line_half],[0,num-tr],':',color=grid_color2,linewidth=lw)
	plt.plot([line,line],[-100,num-tr],'-',color=grid_color1,linewidth=lw)
	plt.plot([-100,num-tr],[line,line],'-',color=grid_color1,linewidth=lw)

for i in xrange(7):
	plt.text(np.sum(num_x_type[0:i+1])+halves[i+1]-tr-25,-50, types[i], fontsize=12,fontweight='bold',color=colors[i])
	plt.text(-80,np.sum(num_x_type[0:i+1])+halves[i+1]-tr+10, types[i], fontsize=12,fontweight='bold',color=colors[i])



