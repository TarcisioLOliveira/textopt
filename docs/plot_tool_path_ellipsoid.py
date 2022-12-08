#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

vc_init = 60
f_uet_init = 20000
Ax_uet_init = 2
Az_uet_init = 5
clear_init = 70

def update_plot_param(vc, f_uet, Ax_uet, Az_uet, clear):
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

def xt(t, vc, f_uet, Ax_uet, Az_uet):
    return -Ax_uet*np.sin(2*math.pi*f_uet*t) + vc*t

def update_plot(vc, f_uet, Ax_uet, Az_uet, clear):
    vc *= 1e6/60

    xoffset_uet = 0
    delta_uet = vc/f_uet

    sq2 = math.sqrt(2)

    tanc = math.tan(clear*math.pi/180)
    v_crit = 2*math.pi*f_uet*Ax_uet
    d_crit = v_crit/f_uet


    X = np.linspace(0, 2*delta_uet, 1000)
    Z = np.linspace(0, 2*delta_uet, 1000)

    if vc <= v_crit:
        t1 = math.acos(vc/v_crit)/(2*math.pi*f_uet)
        t2 = 1.0/f_uet - t1

        # h1 = Az_uet*math.sin((math.pi*math.pi*Ax_uet)/delta_uet) + Az_uet
        # h2 = 2*Az_uet - h1
        # print(h1)

        # Works better with half-senoid
        # H = Az_uet*v_crit/vc
        # h1 =  H + Az_uet 
        # h2 = -H + Az_uet
        # K = delta_uet*(h1-h2)/(2*(h1+h2))
        # delta_1 =  K + delta_uet/2
        # delta_2 = -K + delta_uet/2

        # Works better with half-ellipse
        #h_c = Az_uet*math.cos(-Ax_uet*math.sin(math.acos(vc/v_crit))/vc + math.acos(vc/v_crit)) + Az_uet
        H = Az_uet*math.sin(math.pi*vc/(2*v_crit))
        h1 =  H + Az_uet
        h2 = -H + Az_uet
        K = delta_uet*(h1+h2)/(2*(h1-h2))
        delta_1 =  K + delta_uet/2
        delta_2 = -K + delta_uet/2

        for i in range(len(X)):
            xcirc = X[i] + xoffset_uet
            mult_uet = math.floor(xcirc / delta_uet)
            x_uet = xcirc - (mult_uet + 0.5)*delta_uet

            z_uet = h1*(1.0 - math.sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)))
            #z_uet = h1*(1.0 - math.cos(math.pi*x_uet/delta_1))

            x_uet = xcirc + 0.5*delta_uet - math.floor(xcirc / delta_uet + 0.5)*delta_uet;
            z_clear = tanc*(delta_uet - x_uet);

            Z[i] = min(z_uet, z_clear);
    else:
        H = Az_uet*math.sin(math.pi*v_crit/(2*vc))#Az_uet*v_crit/vc
        h1 =  H + Az_uet
        h2 = -H + Az_uet
        K = delta_uet*(h1-h2)/(2*(h1+h2))
        delta_1 =  K + delta_uet/2
        delta_2 = -K + delta_uet/2
        for i in range(len(X)):
            xcirc = X[i] + xoffset_uet
            mult_uet = math.floor(xcirc / delta_uet)
            x_uet = xcirc - (mult_uet + 0.5)*delta_uet
            x_uet_z = x_uet + delta_uet/2

            if x_uet < -delta_1/2:
                x_uet2 = x_uet + delta_1/2 + delta_2/2

                z_uet1 = h1 + h2*math.cos(math.pi*x_uet2/delta_2)
                z_uet2 = h1 + h2*math.sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2))

            elif x_uet < delta_1/2:
                x_uet2 = x_uet

                z_uet1 = h1*(1.0 - math.cos(math.pi*x_uet2/delta_1))
                z_uet2 = h1*(1.0 - math.sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)))
            else:
                x_uet2 = x_uet - delta_1/2 - delta_2/2

                z_uet1 = h1 + h2*math.cos(math.pi*x_uet2/delta_2)
                z_uet2 = h1 + h2*math.sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2))

            # z_uet = z_uet1*(math.cos(math.pi*v_crit/(2*vc))) + z_uet2*(1 - math.cos(math.pi*v_crit/(2*vc)))
            z_uet = z_uet1*math.sqrt(1 - (v_crit/vc)**2) + z_uet2*(1 - math.sqrt(1 - (v_crit/vc)**2))

            x_uet = xcirc + 0.5*delta_uet - math.floor(xcirc / delta_uet + 0.5)*delta_uet;
            z_clear = tanc*(delta_uet - x_uet);

            Z[i] = min(z_uet, z_clear);

    return X, Z

fig, ax = plt.subplots()
x2, z2 = update_plot_param(vc_init, f_uet_init, Ax_uet_init, Az_uet_init, clear_init)
line2, = plt.plot(x2, z2, 'r-', lw=2)
x, z = update_plot(vc_init, f_uet_init, Ax_uet_init, Az_uet_init, clear_init)
line1, = plt.plot(x, z, 'b-', lw=2)
xmin = np.min(x2)
xmax = np.max(x2)
zmin = np.min(z2)
zmax = np.max(z2)

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
    valmin=0,
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
    x2, z2 = update_plot_param(vc_slider.val, f_uet_slider.val, Ax_uet_slider.val, Az_uet_slider.val, clear_slider.val)
    line2.set_xdata(x2)
    line2.set_ydata(z2)
    x, z = update_plot(vc_slider.val, f_uet_slider.val, Ax_uet_slider.val, Az_uet_slider.val, clear_slider.val)
    line1.set_xdata(x)
    line1.set_ydata(z)
    xmin = np.min(x2)
    xmax = np.max(x2)
    zmin = np.min(z2)
    zmax = np.max(z2)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([zmin, zmax])
    fig.canvas.draw_idle()

vc_slider.on_changed(update)
f_uet_slider.on_changed(update)
Ax_uet_slider.on_changed(update)
Az_uet_slider.on_changed(update)
clear_slider.on_changed(update)

plt.show()
