#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

f_init = 50
ap_init = 50

alpha1_init = 60
alpha2_init = 60

r_init = 10

def update_plot(f, ap, alpha1, alpha2, r):
    alpha1 *= math.pi/180
    alpha2 *= math.pi/180

    tan1 = math.tan(alpha1)
    tan2 = math.tan(alpha2)

    # Line parameters for tool slopes.
    y1 = -math.sqrt((tan1*r*tan1*r)/(tan1*tan1+1))
    y2 =  math.sqrt((tan2*r*tan2*r)/(tan2*tan2+1))

    # Value of z for y = 0 (with circle center at y = 0).
    b1off = -math.sqrt(r*r - y1*y1) + r + tan1*y1
    b2off = -math.sqrt(r*r - y2*y2) + r - tan2*y2

    # Where the slope lines intersect, so that the height at the edges is the
    # same for both
    yc = (tan2*f+b2off-b1off)/(tan1+tan2)

    # Recalculate yc if it gets into the radius, as that skews the centering
    s1 = -yc
    s2 = f-yc
    if alpha1 < alpha2 and s2 < y2:
        # Slope 2 must be disregarded.
        # Height of slope 1 must be equal to height of tool radius at the
        # other edge.
        a = tan1*tan1 + 1
        b = -2*(tan1*(r-b1off)+f)
        c = b1off*(b1off-2*r)+f*f
        yc = (-b - math.sqrt(b*b - 4*a*c))/(2*a)
    elif alpha2 < alpha1 and s1 > y1:
        # Slope 1 must be disregarded.
        # Height of slope 2 must be equal to height of tool radius at the
        # other edge.
        a = tan2*tan2 + 1
        b = -(2*tan2*(b2off-r)+2*tan2*tan2*f)
        c = tan2*tan2*f*f + 2*tan2*f*(b2off-r) + b2off*(b2off-2*r)
        yc = (-b + math.sqrt(b*b - 4*a*c))/(2*a)

    # Calculate whether the angles of the cutting tool intercept the surface
    # first or are stopped by the feed rate.
    #
    # Check the intersection with the feed rate edges. If they're greater than
    # zero, the cut line is narrower than the feed line, so max height is
    # zero. Otherwise, it's the height at the intersection
    #
    # Implementation of this section is based on `texture_shallow`
    h = 0
    s1 = -yc
    s2 = f-yc
    if s1 > y1 and s2 < y2:
        # Width falls entirely within tool radius.
        yc = f/2
        h = -math.sqrt(r*r - yc*yc) + r
    elif s1 > y1:
        h = tan2*(f-yc) + b2off
    elif s2 < y2:
        h = tan1*yc + b1off
    else:
        # Width does not intersect the radius at all.
        h = tan1*yc + b1off

    
    min_z = -ap
    # If it's greater than zero, max_z must be zero
    max_z = min(0.0, h+min_z)

    # Value of z for y = 0 (global)
    line_root1 = -math.sqrt(r*r - y1*y1) + r - ap + tan1*(yc+y1)
    line_root2 = -math.sqrt(r*r - y2*y2) + r - ap - tan2*(yc+y2)

    y = np.linspace(0, f, 1000)

    Y1 = []
    Y2 = []
    Y3 = []
    Z1 = []
    Z2 = []
    Z3 = []

    for i in range(len(y)):
        if y[i] <= y1 + yc:
            Z1.append(-tan1*y[i] + line_root1)
            Y1.append(y[i])
            Z1[-1] = min(Z1[-1], max_z)
        elif y[i] <= y2 + yc:
            yy = y[i] - yc
            Z2.append(-math.sqrt(r*r - yy*yy) + r - ap)
            Y2.append(y[i])
            Z2[-1] = min(Z2[-1], max_z)
        else:
            Z3.append(tan2*y[i] + line_root2)
            Y3.append(y[i])
            Z3[-1] = min(Z3[-1], max_z)

    return Y1, Y2, Y3, Z1, Z2, Z3

fig, ax = plt.subplots()
y1, y2, y3, z1, z2, z3 = update_plot(f_init, ap_init, alpha1_init, alpha2_init, r_init)
line1, = plt.plot(y1, z1, 'b-', lw=2)
line2, = plt.plot(y2, z2, 'r-', lw=2)
line3, = plt.plot(y3, z3, 'b-', lw=2)
ax.set_xlim([np.min(y1), np.max(y3)])
ax.set_ylim([np.min(z2), np.max([z1, z3])])
ax.set_xlabel('Y [um]')
ax.set_ylabel('Z [um]')

plt.subplots_adjust(left=0.50, right=0.97)

axf = plt.axes([0.15, 0.9, 0.2, 0.03])
f_slider = Slider(
    ax=axf,
    label='f [um/rev]',
    valmin=1,
    valmax=200,
    valinit=f_init,
)
axap = plt.axes([0.15, 0.8, 0.2, 0.03])
ap_slider = Slider(
    ax=axap,
    label='ap [um]',
    valmin=1,
    valmax=200,
    valinit=ap_init,
)
axalpha1 = plt.axes([0.15, 0.7, 0.2, 0.03])
alpha1_slider = Slider(
    ax=axalpha1,
    label='alpha1 [deg]',
    valmin=1,
    valmax=89,
    valinit=alpha1_init,
)
axalpha2 = plt.axes([0.15, 0.6, 0.2, 0.03])
alpha2_slider = Slider(
    ax=axalpha2,
    label='alpha2 [deg]',
    valmin=1,
    valmax=89,
    valinit=alpha2_init,
)
axr = plt.axes([0.15, 0.5, 0.2, 0.03])
r_slider = Slider(
    ax=axr,
    label='r [um]',
    valmin=1,
    valmax=200,
    valinit=r_init,
)

def update(val):
    y1, y2, y3, z1, z2, z3 = update_plot(f_slider.val, ap_slider.val, alpha1_slider.val, alpha2_slider.val, r_slider.val)
    line1.set_xdata(y1)
    line1.set_ydata(z1)
    line2.set_xdata(y2)
    line2.set_ydata(z2)
    line3.set_xdata(y3)
    line3.set_ydata(z3)
    if len(y1) > 0:
        ymin = np.min([np.min(y1), np.min(y2)])
        zmax = np.max([np.max(z1), np.max(z2)])
    else:
        ymin = np.min(y2)
        zmax = np.max(z2)
    if len(y3) > 0:
        ymax = np.max([np.max(y2), np.max(y3)])
        zmax = np.max([zmax, np.max(z3)])
    else:
        ymax = np.max(y2)

    ax.set_xlim([ymin, ymax])
    zmin = np.min(z2)
    ax.set_ylim([zmin, zmax])
    fig.canvas.draw_idle()

f_slider.on_changed(update)
ap_slider.on_changed(update)
alpha1_slider.on_changed(update)
alpha2_slider.on_changed(update)
r_slider.on_changed(update)

plt.show()
