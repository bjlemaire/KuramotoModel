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

found_min=0
for count, tthh in enumerate(ths):
    if tthh>0.1 and not found_min:
        found_min=1
        minth = count
    if tthh==1.0:
        maxth = count


listt_s = [0.786568, 0.840939, 0.754871, 0.609112, 0.627880, 0.756211, 0.699419]
listt_s_2 = [2.278479, 2.278221, 2.274414, 2.291461, 2.275544, 2.278048]
t_ss = np.mean(listt_s)
t_min = np.min(listt_s)
t_max = np.max(listt_s)
t_ss_2 = np.mean(listt_s_2)
t_min_2 = np.min(listt_s_2)
t_max_2 = np.max(listt_s_2)
if maxth==100:
    print("All good.")
print("minth=",minth)
print("maxth=",maxth)
minth=51
thl = [np.log(itm) for itm in ths[minth:maxth]]
tsl = [np.log(itm) for itm in ts[minth:maxth]]
A = np.vstack([thl, np.ones(len(thl))]).T
a,b = np.linalg.lstsq(A, tsl, rcond=None)[0]
print("t=alpha*th^beta with alpha={} and beta={}".format(np.exp(b),a))
xx = np.linspace(ths[30], ths[-1],100)
yy = [np.exp(b)*np.power(itm,a) for itm in xx]
fig, ax = plt.subplots(1,figsize=(8,6))
ax.loglog([0.01,2],[t_ss_2,t_ss_2],'-',color='#00358c', alpha=0.9, label = 'Standard Method, 4 cores')
ax.loglog([0.01,2],[t_ss,t_ss],'-',color='#1061e5', alpha=0.9, label = 'Standard Method, 72 cores')
ax.fill_between([0.01,2],[t_max_2,t_max_2],[t_min_2,t_min_2],color='#00358c',alpha=0.3)
ax.fill_between([0.01,2],[t_max,t_max],[t_min,t_min],color='#1061e5',alpha=0.3)
ax.loglog(xx,yy,'-',color = '#05b71d', alpha=1.0, label=r'$\alpha\cdot \theta_{bh}^\beta$ : $\alpha\simeq0.155$, $\beta\simeq -1.04$') # #1532ef
ax.loglog(ths,ts,'-x', color='#d30a0a', label='Datapoints - BH implem., 72 cores')
ax.fill_between(ths, mint, maxt, color = '#d30a0a', alpha=0.3)
plt.legend(loc='lower left', fancybox=True, frameon=True, facecolor='white',fontsize=12)
plt.grid(which='both')
ax.set_xlabel(r"$\theta_{bh}$", fontweight='bold', fontsize=15)
ax.set_ylabel(r"Execution Time", fontweight = 'bold', fontsize=13)
ax.set_xlim([0.01,2])
ax.set_ylim([0.09,4])
ax.axvline(x=0.5,linewidth=1.5,color='#2d2d2d')
ax.text(0.29,0.12,s=r'$\theta_{bh}=0.5$')
ax.set_title(r'Execution Time as a function of $\theta_{Barnes-Hut}$'+ '\n'+r'for $J=0.1$, $K=1$, $T=1.0$, and fixed timestep'+ r' $h=0.05$.', fontweight='bold', fontsize=14)
plt.savefig("ThetaExecTime.png", dpi=600)
plt.show()
