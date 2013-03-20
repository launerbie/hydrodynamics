#!/usr/bin/env python
from amuse.lab import *
from support import read_from_hdf5 as read
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib as mpl 

def main():
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.5e10, 0.5e10), ylim=(-0.5e10, 0.5e10))
    #ax = plt.axes(xlim=(-1.5,1.5), ylim=(-1.5,1.5))
    line, = ax.plot([], [], **yellowdots)

    def genex():
        results = read('N1000n1000.hdf5')
        x = results.positions.x.value_in(units.m)
        y = results.positions.y.value_in(units.m)
        i=0
        length = len(x)
        while i < length:
            yield (x[i], y[i])
            i +=1

    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        line.set_data(i)
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=genex, interval=1, blit=True)
    plt.show()

def rundark():
    mpl.rc('lines', linewidth=1, color='w')
    mpl.rc('patch', edgecolor='w')
    mpl.rc('text',  color='w')
    mpl.rc('font',  size=9, family='sans-serif')
    mpl.rc('axes', facecolor='k', edgecolor='w', labelcolor='w', \
            color_cycle=[ 'w','r','g','y',  'c', 'm', 'b', 'k'],\
            labelsize=9)
    mpl.rc('xtick', color='w')
    mpl.rc('ytick', color='w')
    mpl.rc('grid', color='w')
    mpl.rc('figure', facecolor='k', edgecolor='k')
    mpl.rc('savefig', facecolor='k', edgecolor='k')

def runbright():
    mpl.rc('lines', linewidth=1, color='w')
    mpl.rc('patch', edgecolor='w')
    mpl.rc('text',  color='k')
    mpl.rc('font',  size=9, family='sans-serif')
    mpl.rc('axes', facecolor='w', edgecolor='k', labelcolor='k', \
            color_cycle=[ 'k','r','g','y',  'c', 'm', 'b', 'w'],\
            labelsize=9)
    mpl.rc('xtick', color='k')
    mpl.rc('ytick', color='k')
    mpl.rc('grid', color='k')
    mpl.rc('figure', facecolor='w', edgecolor='w')
    mpl.rc('savefig', facecolor='w', edgecolor='w')

if __name__:
    dark=True
    redline     = dict(c='r', ls="-", lw=1, alpha=1.0)
    yellowline  = dict(c='y', ls="-", lw=1, alpha=1.0)
    blueline    = dict(c='b', ls="-", lw=1, alpha=1.0)
    yellowdots  = dict(c='y', ls="o", mfc="y", mec="y", \
                  marker='o', alpha=1.0, ms=2)
    reddots     = dict(c='r', ls="o", mfc="r", mec="r", \
                  marker='o', alpha=1.0, ms=1)
    greendots   = dict(c='g', ls="o", mfc="g", mec="g", \
                  marker='o', alpha=1.0, ms=1)
    magentadots = dict(c='m', ls="o", mfc="m", mec="m", \
                  marker='o', alpha=1.0, ms=1)
    cyandots    = dict(c='c', ls="o", mfc="c", mec="c", \
                  marker='o', alpha=1.0, ms=1)
    bluedots  = dict(c='b', ls="o", mfc="b", mec="b", \
                  marker='o', alpha=1.0, ms=1)
    if dark == True:
        rundark()
        line = dict(c='w', ls="-", lw=1, alpha=1.0)
        dots = dict(c='w', ls="o", mfc="w", mec="w", \
                      marker='o', alpha=1.0, ms=1)
        errdots = dict(fmt='o', ls="o", ecolor='r', alpha=1.0)
    else:
        runbright()
        line = dict(c='k', ls="-", lw=1, alpha=1.0)
        dots = dict(c='k', ls="o", mfc="k", mec="k", \
                      marker='o', alpha=1.0, ms=1)
        errdots = dict(fmt='o', ls="o", ecolor='r', alpha=1.0)
    main()
