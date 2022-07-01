# genera punti distribuiti NON uniformemente su una superficie sferica.
# mostra che l'estrazione uniforme di theta e phi non è sufficiente
# a garantire una distribuzione uniforme su tutta la sfera

import random as rd
import numpy as np
from matplotlib import pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for i in range(0, 500):
    th=rd.uniform(0, np.pi)
    phi=rd.uniform(0, 2*np.pi)

    xx = np.sin(th)*np.cos(phi)
    yy = np.sin(th)*np.sin(phi)
    zz = np.cos(th)
    ax.scatter(xx,yy,zz, color = [1-(zz+1)/2,0,(zz+1)/2])

plt.xlabel("x")
plt.xlabel("y")
plt.xlabel("z")

plt.show()
