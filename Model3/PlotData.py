import numpy as np
import matplotlib.pyplot as plt

with open("RangeTheta5.dat","r") as temp_file:
    data = temp_file.readlines()

ths = []
ts = []
for line in data:
    th, t1, t2, t3, t4, t5 = [float(itm) for itm in line.split()]
    ths.append(th)
    ts.append(np.mean([t1,t2,t3,t4,t5]))

plt.plot(ths,ts,'-+')
plt.show()
