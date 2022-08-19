#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

file = open("plot_vc.txt", "r")

vc = file.readline().split(" ");
vci = float(vc[0])
vcf = float(vc[1])
vcd = float(vc[2])
vc = np.linspace(vci, vcf, int((vcf - vci)/vcd)+1)

A = np.zeros([len(vc)])
Sa = np.zeros([len(vc)])
for vci in range(len(vc)):
    points = file.readline().split(" ")
    A[vci] = float(points[0])
    Sa[vci] = float(points[1])

file.close()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlabel("vc [m/min]")
ax1.set_ylabel("A [um^2]")
ax1.plot(vc, A)
           

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel("vc [m/min]")
ax2.set_ylabel("Sa [um]")
ax2.plot(vc, Sa)

plt.show()
