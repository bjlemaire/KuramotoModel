import numpy as np
import matplotlib.pyplot as plt

with open("RangeTheta.dat","r") as temp_file:
    data = temp_file.readlines()

ths = []
ts = []
mint = []
maxt = []
for line in data:
    th, t1, t2, t3, t4, t5 = [float(itm) for itm in line.split()]
    ths.append(th)
    ts.append(np.mean([t1,t2,t3,t4,t5]))
    mint.append(np.min([t1,t2,t3,t4,t5]))
    maxt.append(np.max([t1,t2,t3,t4,t5]))

thl = [np.log(itm) for itm in ths[59:]]
tsl = [np.log(itm) for itm in ts[59:]]
A = np.vstack([thl, np.ones(len(thl))]).T
a,b = np.linalg.lstsq(A, tsl, rcond=None)[0]
print("t=alpha*th^beta with alpha={} and beta={}".format(np.exp(b),a))
xx = np.linspace(ths[30], ths[-1],100)
yy = [np.exp(b)*np.power(itm,a) for itm in xx]
fig, ax = plt.subplots(1,figsize=(8,6))
ax.loglog(xx,yy,'-',color = '#05b71d', alpha=1.0, label=r'$\alpha\cdot \theta_{bh}^\beta$ with $\alpha\simeq0.19$, $\beta\simeq -1.2$') # #1532ef
ax.loglog(ths,ts,'-x', color='#d30a0a', label='Datapoints')
ax.fill_between(ths, mint, maxt, color = '#d30a0a', alpha=0.3)
plt.legend(loc='lower left', fancybox=True, frameon=True, facecolor='white',fontsize=12)
plt.grid(which='both')
ax.set_xlabel(r"$\theta_{bh}$", fontweight='bold', fontsize=15)
ax.set_ylabel(r"Execution Time", fontweight = 'bold', fontsize=13)
ax.set_xlim([0.01,0.5])
ax.set_ylim([0.4,8])
ax.set_title(r'Execution Time as a function of $\theta_{Barnes-Hut}$'+ '\n'+r'for $J=0.1$, $K=1$, $T=1.0$, and fixed timestep'+ r' $h=0.05$.', fontweight='bold', fontsize=14)
plt.savefig("ThetaExecTime.png", dpi=600)
plt.show()
