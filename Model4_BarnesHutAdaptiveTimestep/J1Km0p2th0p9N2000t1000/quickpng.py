import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FuncFormatter, MultipleLocator

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
        plt.plot(x,y, '.', color = rgb(θ))
        listx.append(x)
        listy.append(y)
        listth.append(np.fmod(θ, 2*np.pi))
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
    plt.savefig(file[:-4]+'.png', dpi=300)
    plt.close('all')
    fig, ax = plt.subplots(1, figsize = (5.5,4.5))
    plt.plot(listphi,listth, '.', color = '#6603a0', alpha =0.8)
    plt.xlabel(r'$\phi$', fontsize = 13)
    plt.ylabel(r'$\theta$', fontsize = 13, rotation=0)
    plt.xlim([0, 2*np.pi])
    plt.ylim([0, 2*np.pi])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
    ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
    ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi))
    ticklab = ax.xaxis.get_ticklabels()[0]
    trans = ticklab.get_transform()
    ax.xaxis.set_label_coords(6.5, 0.05, transform=trans)
    ax.yaxis.set_label_coords(-0.01,1.02)
    plt.title('t={:05.3f}'.format(t), fontweight = 'bold', fontsize = 14)
    #plt.tight_layout()
    plt.savefig(file[:-4]+'_phase.png',dpi=300)
    plt.close('all')
    listWp.append(Wp)
    listWm.append(Wm)

fig, ax = plt.subplots(1, figsize = (6.5,4.5))
plt.plot(listT,listWp,'.', color = '#09a02c', label = r'$W_p$')
plt.plot(listT,listWm,'.', color = '#a01f09', label = r'$W_m$')
plt.legend(fancybox=True, frameon = True, facecolor = 'white')
plt.xlabel("Time", fontsize = 13)
plt.savefig('WpWm.png',dpi=300)
plt.close('all')
    
