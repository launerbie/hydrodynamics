#!/usr/bin/env python

import matplotlib as mpl

def rundark():
    mpl.rc('lines', linewidth=1, color='w')
    mpl.rc('patch', edgecolor='w')
    mpl.rc('text', color='w')
    mpl.rc('font', size=9, family='sans-serif')
    mpl.rc('axes', facecolor='k', edgecolor='w', labelcolor='w',\
            color_cycle=[ 'w','r','g','y', 'c', 'm', 'b', 'k'],\
            labelsize=9)
    mpl.rc('xtick', color='w')
    mpl.rc('ytick', color='w')
    mpl.rc('grid', color='w')
    mpl.rc('figure', facecolor='k', edgecolor='k')
    mpl.rc('savefig', dpi=150, facecolor='k', edgecolor='k')

def runbright():
    mpl.rc('lines', linewidth=1, color='w')
    mpl.rc('patch', edgecolor='w')
    mpl.rc('text', color='k')
    mpl.rc('font', size=9, family='sans-serif')
    mpl.rc('axes', facecolor='w', edgecolor='k', labelcolor='k', \
            color_cycle=[ 'k','r','g','y', 'c', 'm', 'b', 'w'],\
            labelsize=9)
    mpl.rc('xtick', color='k')
    mpl.rc('ytick', color='k')
    mpl.rc('grid', color='k')
    mpl.rc('figure', facecolor='w', edgecolor='w')
    mpl.rc('savefig', dpi=150, facecolor='w', edgecolor='w')

if __name__ == "theme":
    dark=True
    redline     = dict(c='r', ls="-", lw=1, alpha=1.0)
    yellowline  = dict(c='y', ls="-", lw=1, alpha=1.0)
    blueline    = dict(c='b', ls="-", lw=1, alpha=1.0)
    magentaline = dict(c='m', ls="-", lw=1, alpha=1.0)
    cyanline    = dict(c='c', ls="-", lw=1, alpha=1.0)
    greenline   = dict(c='g', ls="-", lw=1, alpha=1.0)
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

