import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

f = lambda θ: 0.45*(1+np.cos(θ))
#f = lambda θ: 0.25*(1+np.cos(θ))

def rgb(θ):
#    return [f(θ-(2/3)), 0.15, f(1.2*θ+(2/3))] 
    return [f(θ), f(θ-(3/3)), f(θ+(1/3))] 
#    return [f(θ), f(θ-(2/3)), f(θ+(1/3))] 
#    return [0.2, 0.4, f(θ+(1/3))] 

dense_files = []
for file in list(os.listdir('.')):
    if file.startswith('Dense') and file.endswith('.txt'):
        dense_files.append(file)
dense_files = sorted(dense_files)

listT = []
listWp = []
listWm = []

myTan = lambda y,x: np.arctan2(y,x) if np.arctan2(y,x)>0 else 2*np.pi+np.arctan2(y,x)

for count, file in enumerate(dense_files):
    fig, ax = plt.subplots(1, figsize = (5,5))
    listx = []
    listy = []
    listth = []
    listphi = []
    Wp = 0.0
    Wm = 0.0
    with open(file, 'r') as temp_file:
        data = temp_file.readlines()
    listT.append(float(data[0].split()[0]))
    for line in data:
        t, x, y, θ = [float(itm) for itm in line.split()]
        ø = (myTan(y,x))
        plt.plot(x,y, 'o', color = rgb(θ))
        listx.append(x)
        listy.append(y)
        listth.append(θ)
        listphi.append(ø)
        Wp += (1/len(data))*np.exp(complex(0,(ø+θ)))
        Wm += (1/len(data))*np.exp(complex(0,(ø-θ)))
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
    fig, ax = plt.subplots(1, figsize = (5,5))
    plt.plot(listphi,listth, '.', color = '#2a5fc1')
    plt.xlabel(r'$\phi$', fontsize = 13)
    plt.ylabel(r'$\theta$', fontsize = 13)
    plt.xlim([0, 2*np.pi])
    plt.ylim([0, 2*np.pi])
    plt.title('t={:05.3f}'.format(t), fontweight = 'bold', fontsize = 14)
    plt.savefig(file[:-4]+'_phase.png')
    plt.close('all')
    listWp.append(Wp)
    listWm.append(Wm)

fig, ax = plt.subplots(1, figsize = (5,5))
plt.plot(listT,listWp,'.', color = '#09a02c', label = r'$W_p$')
plt.plot(listT,listWm,'.', color = '#a01f09', label = r'$W_m$')
plt.legend(loc = 'upper left', fancybox=True, frameon = True, facecolor = 'white')
plt.savefig('WpWm.png')
plt.close('all')
    
