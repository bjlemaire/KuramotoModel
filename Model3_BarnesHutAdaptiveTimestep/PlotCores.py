import matplotlib.pyplot as plt
import numpy as np

with open("RangeCores.dat","r") as temp_file:
    data = temp_file.readlines();

ncore = []
ts = []
mint = []
maxt = []
for line in data:
    ncore.append(int(line.split()[0]))
    ts.append(np.mean([float(itm) for itm in line.split()[1:]]))
    mint.append(np.min([float(itm) for itm in line.split()[1:]]))
    maxt.append(np.max([float(itm) for itm in line.split()[1:]]))
spd = [ts[-1]/itm for itm in ts[:-1]]
minspd = [mint[-1]/itm for itm in maxt[:-1]]
maxspd = [maxt[-1]/itm for itm in mint[:-1]]
fig, ax = plt.subplots(1, figsize=(7,5))
ax.plot(ncore, ts,'-x',color='#af0f3a')
ax.fill_between(ncore, mint, maxt, color='#af0f3a', alpha=0.3)
ax.set_xlabel(r'Number of Cores')
ax.set_ylabel(r'Execution Time (s)')
ax.set_title('Execution Time as a function of the number of computing cores,\n for'+r' $J=0.1$, $K=1$, $\theta=0.5$ on a c5.18large AWS instance.')
plt.grid(which='both')
ax.set_xlim([0,40])
ax.set_ylim([3,40])
plt.savefig("ExecTimeCore.png",dpi=300)
fig, ax = plt.subplots(1, figsize=(7,5))
ax.plot(ncore[:-1], spd,'-x',color='#25af0f')
ax.fill_between(ncore[:-1], minspd, maxspd, color='#25af0f', alpha=0.3)
ax.set_xlabel(r'Number of Cores')
ax.set_ylabel(r'Speedup')
ax.set_title('Speedup as a function of the number of computing cores,\nfor'+r' $J=0.1$, $K=1$, $\theta=0.5$ on a c5.18large AWS instance.')
ax.set_xlim([2, 36])
plt.grid(which='both')
plt.savefig("SpeedupCore.png",dpi=300)
plt.show()
