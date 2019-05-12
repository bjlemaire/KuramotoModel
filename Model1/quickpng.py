import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

f = lambda θ: 0.45*(1+np.cos(θ))

def rgb(θ):
    return [f(θ), f(θ-(2/3)), f(θ+(2/3))] 

dense_files = []
for file in list(os.listdir('.')):
    if file.startswith('DenseQ') and file.endswith('_iii.txt'):
        dense_files.append(file)
dense_files = sorted(dense_files)

for count, file in enumerate(dense_files):
    fig, ax = plt.subplots(1, figsize = (5,5))
    listx = []
    listy = []
    with open(file, 'r') as temp_file:
        data = temp_file.readlines()
    for line in data:
        t, x, y, θ = [float(itm) for itm in line.split()]
        plt.plot(x,y, 'o', color = rgb(θ))
        listx.append(x)
        listy.append(y)
    minx = min(listx)
    miny = min(listy)
    maxx = max(listx)
    maxy = max(listy)
    lenx = (maxx-minx)
    leny = (maxy-miny)
    lim_min = min([minx-0.05*lenx, miny-0.05*leny])
    lim_max = max([maxx+0.05*lenx, maxy+0.05*leny])
    plt.xlim([lim_min, lim_max])
    plt.ylim([lim_min, lim_max])
    plt.title('t={:05.3f}'.format(t), fontweight = 'bold', fontsize = 14)
    plt.savefig(file[:-4]+'.png')
    plt.close('all')
