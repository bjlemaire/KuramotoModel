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


plt.figure()
plt.loglog(NN, ts)
plt.loglog(NN, tbh)
plt.grid(which='both')
plt.show()
