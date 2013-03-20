#!/usr/bin/env python
import sys
import h5py

from optparse import OptionParser

from amuse.units import units

from matplotlib import pyplot as plt
import matplotlib as mpl

from support import read_from_hdf5

def main(options):
    filename = options.hdf5file
    results = read_from_hdf5(filename)
    return 0 

def plot_Ndependence(results):
    pass

def plot_steptimes(steptimes):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    ax1.plot(range(len(steptimes)), steptimes, **yellowline)
    ax1.set_xlabel('N')
    ax1.set_ylabel('step time [in seconds]')
    plt.show()

def plot_energy(energy, energy_error, value_in="units.J"):
    energy = eval("energy.value_in("+value_in+")")
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.plot(range(len(energy)), energy, **yellowline)
    ax1.set_xlabel('N')
    ax1.set_ylabel('energy [%s]'%value_in)

    ax2.plot(range(len(energy)), energy_error, **yellowline)
    ax2.set_xlabel('N')
    ax2.set_ylabel('energy_error')
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

if __name__ == "plotter":
    dark=True
    redline     = dict(c='r', ls="-", lw=1, alpha=1.0)
    yellowline  = dict(c='y', ls="-", lw=1, alpha=1.0)
    blueline    = dict(c='b', ls="-", lw=1, alpha=1.0)
    yellowdots  = dict(c='y', ls="o", mfc="y", mec="y", \
                  marker='o', alpha=1.0, ms=1)
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

def parse_sysargs(sysargs):
    """ Yeah, it's an OptionParser. """
    parser = OptionParser()
    parser.set_defaults(test=False, skip_table=False)


    parser.add_option("-f", "--file", action="store",\
                      dest="hdf5file", default=None,\
                      help="HDF5 file to open")

    parser.add_option("-e", "--energy", action="store_true",\
                      dest="plot_NvsE", default=False,\
                      help="help txt")

    parser.add_option("-t", "--times", action="store_true", \
                      dest="plot_NvsT", default=False,\
                      help="help txt")

    options, args = parser.parse_args(sysargs[1:])
    return options

if __name__ in ('__main__'):
    options = parse_sysargs(sys.argv) 
    main(options)
