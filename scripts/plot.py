#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys

path = ""
if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    print("Missing input file!")
    exit()

file = open(path, "r")

dim = int(file.readline())
if dim == 2:
    label_x = file.readline()

    f = file.readline().split(" ");
    fi = float(f[0])
    ff = float(f[1])
    fd = float(f[2])
    f = np.linspace(fi, ff, int((ff - fi)/fd)+1)

    label_y = file.readline()

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
    ax1.set_xlabel(label_x)
    ax1.set_ylabel(label_y)
    ax1.set_zlabel("A [um^2/um^2]")
    ax1.plot_surface(f, ap, A, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.set_xlabel(label_x)
    ax2.set_ylabel(label_y)
    ax2.set_zlabel("Sa [um]")
    ax2.plot_surface(f, ap, Sa, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    plt.show()

elif dim == 1:

    label_x = file.readline()

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
    ax1.set_xlabel(label_x)
    ax1.set_ylabel("A [um^2/um^2]")
    ax1.plot(vc, A)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel(label_x)
    ax2.set_ylabel("Sa [um]")
    ax2.plot(vc, Sa)

    plt.show()
