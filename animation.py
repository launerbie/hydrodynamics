#!/usr/bin/env python
import sys
from amuse.lab import *
from support import read_from_hdf5 as read
import matplotlib.pyplot as plt
import matplotlib as mpl 
import numpy as np
from optparse import OptionParser

def main(options):
    """ Creates images of the particle positions at each timestep."""
    results = read(options.hdf5file)
    x = results.positions.x.value_in(units.RSun)
    y = results.positions.y.value_in(units.RSun)
    
    totalsteps = len(x)
    nr_particles = len(x[0])
    times = results.time
    end_time = times[-1].value_in(units.day)

    for i in range(totalsteps):
        print i
        fig = plt.figure(figsize=(8,8), dpi=300)
        ax1 = fig.add_subplot(111)

        ax1.plot(x[i],y[i], **yellowdots)
        ax1.set_title("Particles:%i,  Steps:%i,  End-time:%f days,  t:%f days"%(nr_particles, totalsteps,\
                           end_time, times[i].value_in(units.day)))
        ax1.set_xlim(-2.0, 2.0)
        ax1.set_ylim(-2.0, 2.0)
        ax1.set_xlabel('RSun')
        ax1.set_ylabel('RSun')
        
        zerostring = zeroes(i)
        plt.savefig("images/"+zerostring+str(i)+".png", dpi=120)
        fig.clf()
        plt.close()

    return 0

def zeroes(framenr):
    """ Checks how many zeroes need to be added. You don't want 1.jpg 
    but 0001.jpg """
    zeroes=0
    if framenr/1000 == 0:
        zeroes += 1
        if framenr/100 == 0:
            zeroes += 1
            if framenr/10 == 0:
                zeroes += 1
    zerostring = zeroes*"0"
    return zerostring

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

def parse_sysargs(sysargs):
    """ OptionParser. """
    parser = OptionParser()
    parser.set_defaults(test=False, skip_table=False)

    parser.add_option("-f", "--file", action="store",\
                      dest="hdf5file", default=None,\
                      help="HDF5 file to open")

    options, args = parser.parse_args(sysargs[1:])
    return options


if __name__ in ('__main__'):
    dark=True
    redline     = dict(c='r', ls="-", lw=1, alpha=1.0)
    yellowline  = dict(c='y', ls="-", lw=1, alpha=1.0)
    blueline    = dict(c='b', ls="-", lw=1, alpha=1.0)
    yellowdots  = dict(c='y', ls="o", mfc="y", mec="r", \
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


    options = parse_sysargs(sys.argv)
    main(options)

