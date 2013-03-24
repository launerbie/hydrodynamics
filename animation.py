#!/usr/bin/env python
import os
import sys
import shutil
import subprocess
from progressbar import progressbar as pb
from optparse import OptionParser

import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt

from amuse.units import units
from support import read_from_hdf5 as read

def main(options):
    """ Creates animation of the particles by making png's and encode
    the png's using ffmpeg."""
    results = read(options.hdf5file)

    if os.path.exists('images'):
        shutil.rmtree('images')
        os.makedirs('images')
    else:
        os.makedirs('images')

    create_images(results, axis_range=options.axis_range)

    if options.outputfile:
        filename = options.outputfile
    else:
        filename = options.hdf5file[:-5]+"r"+str(options.axis_range)+".mp4"

    create_movie(outputfile=filename, fps=options.fps)
    return 0

def create_movie(outputfile="movie.mp4", fps=30):
    """ Encode a list of png's located in hydrodynamics/images into
    an animation. """
    command = 'ffmpeg -q:v 5 -r %i -b:v 9600 -y -i images/%%04d.png %s'%(\
              fps, outputfile)
    subprocess.call(command, shell=True)
    print "File written to: %s"%outputfile

def create_images(results, axis_range=4.0, scalefactor=40.0):
    """ Batch creates plots (.png) of the particle's positions and
    write them to images."""
    x = results.positions.x.value_in(units.RSun)
    y = results.positions.y.value_in(units.RSun)
    z = results.positions.z.value_in(units.RSun)
    
    totalsteps = len(x)
    nr_particles = len(x[0])
    times = results.time
    end_time = times[-1].value_in(units.day)

    widget = drawwidget("Generating images")
    pbar = pb.ProgressBar(widgets=widget, maxval=totalsteps).start()
    for i in range(totalsteps):
        fig = plt.figure(figsize=(8,8), dpi=300)
        ax1 = fig.add_subplot(111)
  
        size = (scalefactor*z[i])/axis_range**2

        ax1.scatter(x[i],y[i], s=size , edgecolor='b', facecolor='w')
        ax1.set_title("Particles:%i  Steps:%i  End-time:%f days  t:%f days"%(nr_particles, totalsteps,\
                           end_time, times[i].value_in(units.day)))
        ax1.set_xlim((-1)*axis_range, axis_range)
        ax1.set_ylim((-1)*axis_range, axis_range)
        ax1.set_xlabel('RSun')
        ax1.set_ylabel('RSun')
        
        zerostring = zeroes(i)
        plt.savefig("images/"+zerostring+str(i)+".png", dpi=120)
        fig.clf()
        plt.close()
        pbar.update(i)
    pbar.finish()

def drawwidget(proces_discription):
    """ Formats the progressbar. """
    widgets = [proces_discription+": ", pb.Percentage(), ' ',
               pb.Bar(marker='#',left='[',right=']'),
               ' ', pb.ETA()]
    return widgets

def zeroes(framenr):
    """ Checks how many zeroes need to be added. You don't want 1.jpg
    but 0001.jpg """
    #To do: specify how many digits you want. Now it's always 4.
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
    """ Sets the Matplotlib rc parameters to run with a dark
    background."""
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
    """ Sets the Matplotlib rc parameters to run with a bright
    background."""
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
    parser.set_defaults(test=False)

    parser.add_option("-f", "--file", action="store",\
                      dest="hdf5file", default=None,\
                      help="HDF5 file to open")
    parser.add_option("-o", "--output", action="store",\
                      dest="outputfile", default=None,\
                      type='string',\
                      help="Filename for outputted mp4 file.")
    parser.add_option("-F", "--fps", action="store", type='int',\
                      dest="fps", default=30,\
                      help="FPS: frames per second for encoding.")
    parser.add_option("-r", "--range", action="store", type='float',\
                      dest="axis_range", default=4.0,\
                      help="Size of axis in RSun.")
    parser.add_option("-s", "--scale", action="store", type='float',\
                      dest="scalefactor", default=40.0,\
                      help="Scale factor for dots.")

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

