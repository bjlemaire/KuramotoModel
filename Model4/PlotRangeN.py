import matplotlib.pyplot as plt
import numpy as np

with open("RangeN.dat","r") as temp_file:
    data = temp_file.readlines()

NN = []
ts = []
tbh = []

for line in data:
    NN.append(int(line.split()[0]))
    ts.append(float(line.split()[1]))
    tbh.append(float(line.split()[2]))
NNl = np.array([np.log(itm) for itm in NN])
tsl = np.array([np.log(itm) for itm in ts])
tbhl = np.array([np.log(itm) for itm in tbh])
A = np.vstack([NNl, np.ones(len(NNl))]).T
a_s, b_s = np.linalg.lstsq(A, tsl, rcond=None)[0]
a_b, b_b = np.linalg.lstsq(A, tbhl, rcond=None)[0]
print("For standard: alpha={:.8f} and beta={:.6f}".format(np.exp(b_s),a_s))
print("For Barnes-Hut: alpha={:.8f} and beta={:.6f}".format(np.exp(b_b),a_b))
xx= np.linspace(1000,15848,200)

yy_s = [np.exp(b_s)*np.power(itm, a_s) for itm in xx]
yy_b = [np.exp(b_b)*np.power(itm, a_b) for itm in xx]

fig, ax = plt.subplots(1, figsize=(7,5))
ax.loglog(NN, ts, 'x', color = '#2f20bc', alpha = 1.0, label = r'Initial Method')
ax.loglog(xx, yy_s, '-', color = '#2f20bc', alpha=0.7, label=r'$\alpha N^\beta$ - $\alpha\simeq 1.0e-6$, $\beta\simeq 1.823$')
ax.loglog(NN, tbh, 'x', color = '#d66f02', alpha =1.0, label = r'Barnes-Hut Method')
ax.loglog(xx, yy_b, '-', color = '#d66f02', alpha=0.7, label=r'$\alpha N^\beta$ - $\alpha\simeq 5.2e-3$, $\beta\simeq 0.714$')
ax.set_xlabel(r'Number of particles $N$', fontsize = 13)
ax.set_ylabel(r'Execution Time (s)', fontsize = 13)
ax.set_xlim([1000,15868])
ax.set_title("Comparison initial method and Barnes-Hut\nimplementation "+r'with $J=0.1$, $K=1$, $\theta=0.5$', fontsize=14, fontweight='bold')
plt.legend(loc='upper left', frameon=True, fancybox=True, facecolor='w', fontsize=11)
plt.grid(which='both')
plt.savefig("ExecTimeN.png",dpi=300)
plt.show()
