#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys

path = "plot_f_shallow.txt"
if len(sys.argv) > 1:
    path = sys.argv[1]

file = open(path, "r")

f = file.readline().split(" ");
fi = float(f[0])
ff = float(f[1])
fd = float(f[2])
f = np.linspace(fi, ff, int((ff - fi)/fd)+1)

A = np.zeros([len(f)])
Sa = np.zeros([len(f)])
for fi in range(len(f)):
    points = file.readline().split(" ")
    A[fi] = float(points[0])
    Sa[fi] = float(points[1])

file.close()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlabel("f [um/rev]")
ax1.set_ylabel("A [um^2]")
ax1.plot(f, A)
           

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel("f [um/rev]")
ax2.set_ylabel("Sa [um]")
ax2.plot(f, Sa)

plt.show()
