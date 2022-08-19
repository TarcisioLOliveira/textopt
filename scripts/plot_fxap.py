#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys

path = "plot_fxap.txt"
if len(sys.argv) > 1:
    path = sys.argv[1]

file = open(path, "r")

f = file.readline().split(" ");
fi = float(f[0])
ff = float(f[1])
fd = float(f[2])
f = np.linspace(fi, ff, int((ff - fi)/fd)+1)

ap = file.readline().split(" ");
api = float(ap[0])
apf = float(ap[1])
apd = float(ap[2])
ap = np.linspace(api, apf, int((apf - api)/apd)+1)

A = np.zeros([len(f), len(ap)])
Sa = np.zeros([len(f), len(ap)])
for fi in range(len(f)):
    for api in range(len(ap)):
        points = file.readline().split(" ")
        A[fi, api] = float(points[0])
        Sa[fi, api] = float(points[1])

file.close()

f, ap = np.meshgrid(f, ap)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.set_xlabel("f [um/rev]")
ax1.set_ylabel("ap [um]")
ax1.set_zlabel("A [um^2]")
ax1.plot_surface(f, ap, A, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.set_xlabel("f [um/rev]")
ax2.set_ylabel("ap [um]")
ax2.set_zlabel("Sa [um]")
ax2.plot_surface(f, ap, Sa, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()
