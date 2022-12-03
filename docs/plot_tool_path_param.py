#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

vc_init = 60
f_uet_init = 20000
Ax_uet_init = 2
Az_uet_init = 5
clear_init = 10

def update_plot(vc, f_uet, Ax_uet, Az_uet, clear):
    vc *= 1e6/60

    xoffset_uet = 0
    delta_uet = vc/f_uet

    T = 1/f_uet;

    tanc = math.tan(clear*math.pi/180)
    v_crit = 2*math.pi*f_uet*Ax_uet

    t = np.linspace(0, 2*T, 1000)
    X = np.linspace(0, 2*delta_uet, 1000)
    Z = np.linspace(0, 2*delta_uet, 1000)

    X = -Ax_uet*np.sin(2*math.pi*f_uet*t) + vc*t
    Z = Az_uet*np.cos(2*math.pi*f_uet*t) + Az_uet

    return X, Z

fig, ax = plt.subplots()
x, z = update_plot(vc_init, f_uet_init, Ax_uet_init, Az_uet_init, clear_init)
line1, = plt.plot(x, z, 'b-', lw=2)
xmin = np.min(x)
xmax = np.max(x)
zmin = np.min(z)
zmax = np.max(z)

ax.set_xlim([xmin, xmax])
ax.set_ylim([zmin, zmax])
# ax.set_xlim([np.min(y1), np.max(y3)])
# ax.set_xlim([np.min(y1), np.max(y3)])
#ax.set_ylim([np.min(z), np.max([z, z])])
#ax.set_ylim([np.min(z2), np.max([z1, z3])])
ax.set_xlabel('X [um]')
ax.set_ylabel('Z [um]')

plt.subplots_adjust(left=0.50, right=0.97)

axvc = plt.axes([0.15, 0.9, 0.2, 0.03])
vc_slider = Slider(
    ax=axvc,
    label='vc [m/min]',
    valmin=1,
    valmax=200,
    valinit=vc_init,
)
axf_uet = plt.axes([0.15, 0.8, 0.2, 0.03])
f_uet_slider = Slider(
    ax=axf_uet,
    label='f_uet [um]',
    valmin=10000,
    valmax=50000,
    valinit=f_uet_init,
)
axAx_uet = plt.axes([0.15, 0.7, 0.2, 0.03])
Ax_uet_slider = Slider(
    ax=axAx_uet,
    label='Ax_uet [um]',
    valmin=1,
    valmax=100,
    valinit=Ax_uet_init,
)
axAz_uet = plt.axes([0.15, 0.6, 0.2, 0.03])
Az_uet_slider = Slider(
    ax=axAz_uet,
    label='Az_uet [um]',
    valmin=1,
    valmax=100,
    valinit=Az_uet_init,
)
axclear = plt.axes([0.15, 0.5, 0.2, 0.03])
clear_slider = Slider(
    ax=axclear,
    label='clear [um]',
    valmin=1,
    valmax=89,
    valinit=clear_init,
)

def update(val):
    x, z = update_plot(vc_slider.val, f_uet_slider.val, Ax_uet_slider.val, Az_uet_slider.val, clear_slider.val)
    line1.set_xdata(x)
    line1.set_ydata(z)
    xmin = np.min(x)
    xmax = np.max(x)
    zmin = np.min(z)
    zmax = np.max(z)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([zmin, zmax])
    fig.canvas.draw_idle()

vc_slider.on_changed(update)
f_uet_slider.on_changed(update)
Ax_uet_slider.on_changed(update)
Az_uet_slider.on_changed(update)
clear_slider.on_changed(update)

plt.show()
