from time import pthread_getcpuclockid
from turtle import shape
from wsgiref.util import shift_path_info
import matplotlib
import matplotlib.pyplot as plt
from requests import PreparedRequest
import numpy as np
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
import imageio
import sys

C=0

print("Generating image number:")

gen,best = np.loadtxt("bestLen/BL.out", usecols=(0,1), delimiter=' ', unpack='true')

X,Y = np.loadtxt("cit/1citycoord.out",  usecols=(0,1), delimiter=' ', unpack='true')
Best = 500
print("Generating image number (of"+str(Best)+"):")

fig=plt.figure(figsize=(5, 5))
filenames=[]
N=(X.size)

for j in range(0,Best):

    X,Y = np.loadtxt("cit/"+str(j)+"citycoord.out",  usecols=(0,1), delimiter=' ', unpack='true')
    fig=plt.figure(figsize=(8, 8))

    plt.plot(X,Y, label = "Generation "+str(j)+", Lenght ="+str(best[j]))  
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.ylim([-1.5, 1.5])
    plt.xlim([-1.5, 1.5])
    plt.grid(True)
    plt.legend()
    if j%10 == 0:
        print(j, end = ", ")
    plt.savefig("imgs/%s.png"%j,format="png", dpi=150)
    filenames.append("imgs/%s.png"%j)
    plt.close(fig)

with imageio.get_writer('movie.gif', mode='I') as writer:
    i = 0
    print("\nGenerating GIF: ")
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
        i = i+1
        if i%10 == 0:
            print(i, end = ", ")

print("\nEND")