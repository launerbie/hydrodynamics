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
    plot_all(results)
    return 0 

def plot_all(results):
    n_steps = len(results.time)
    N = len(results.positions.x[0])

    lr1 = results.lagrangianradii[:,0]
    lr2 = results.lagrangianradii[:,1]
    lr3 = results.lagrangianradii[:,2]
    lr4 = results.lagrangianradii[:,3]

    Lx = results.angular_momentum.x
    Ly = results.angular_momentum.y
    Lz = results.angular_momentum.z
    L = results.angular_momentum.lengths()

    Ekin = results.kinetic_energy
    Epot = results.potential_energy
    Etot = results.total_energy

    radius_initial = results.radius_initial.as_quantity_in(units.RSun)
    radius_final = results.radius_final.as_quantity_in(units.RSun)

    densities_initial = results.densities_initial
    densities_final = results.densities_final

    times = results.time
   
    fig = plt.figure()

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot( times.value_in(times.unit),     lr1.value_in(lr1.unit), label='10%',  **line)
    ax1.plot( times.value_in(times.unit),     lr2.value_in(lr1.unit), label='25%', **yellowline)
    ax1.plot( times.value_in(times.unit),     lr3.value_in(lr1.unit), label='50%', **redline)
    ax1.plot( times.value_in(times.unit),     lr4.value_in(lr1.unit), label='75%', **greenline)
    ax1.set_xlabel('Time in %s'%times.unit.__str__())
    ax1.set_ylabel('Lagrangian Radius in %s'%lr1.unit.__str__())
    ax1.legend(loc='best')

    ax2.plot( times.value_in(times.unit),    L.value_in(Lx.unit),  label='L',  **line )
    ax2.plot( times.value_in(times.unit),    Lx.value_in(Lx.unit), label='Lx',  **yellowline )
    ax2.plot( times.value_in(times.unit),    Ly.value_in(Lx.unit), label='Ly',  **redline )
    ax2.plot( times.value_in(times.unit),    Lz.value_in(Lx.unit), label='Lz',  **greenline)
    ax2.set_xlabel('Time in %s'%times.unit.__str__())
    ax2.set_ylabel('Angular Momentum in %s'%L.unit.__str__())
    ax2.legend(loc='best')

    ax3.plot(  times.value_in(times.unit),   Ekin.value_in(Ekin.unit),label='Kinetic',  **line)
    ax3.plot(  times.value_in(times.unit),   Epot.value_in(Ekin.unit),label='Potential',  **yellowline)
    ax3.plot(  times.value_in(times.unit),   Etot.value_in(Ekin.unit),label='Total',  **redline)
    #ax3.plot(Eerr.value_in(Ekin.unit), times.value_in(times.unit))    
    ax3.set_xlabel('Time in %s'%times.unit.__str__())
    ax3.set_ylabel('Energy in %s'%Etot.unit.__str__())
    ax3.legend(loc='best')

    ax4.plot(radius_initial.value_in(radius_initial.unit), densities_initial.value_in(densities_initial.unit), label='Initial', **line )
    ax4.plot(radius_final.value_in(radius_final.unit), densities_final.value_in(densities_initial.unit), label='Final', **yellowline )
    ax4.set_xlabel('Radius in %s'%radius_initial.unit.__str__())
    ax4.set_ylabel('Density in %s'%densities_initial.unit.__str__())
    ax4.set_xlim(-2,3)
    ax4.legend(loc='best')
    plt.suptitle('Particles: %i  Steps: %i'%(N, n_steps))
    plt.show()
    pass

#def plot_Ndependence(results):
#    fig = plt.figure()
#    ax1 = fig.add_subplot(211)
#    ax2 = fig.add_subplot(212)
#    
#    ax1.plot(range(len(steptimes)), steptimes, **yellowline)
#    ax1.set_xlabel('N')
#    ax1.set_ylabel('step time [in seconds]')
#    plt.show()
#    
#    fig
#    pass
#
#def plot_steptimes(steptimes):
#    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
#    
#    ax1.plot(range(len(steptimes)), steptimes, **yellowline)
#    ax1.set_xlabel('N')
#    ax1.set_ylabel('step time [in seconds]')
#    plt.show()
#
#def plot_energy(energy, energy_error, value_in="units.J"):
#    energy = eval("energy.value_in("+value_in+")")
#    fig = plt.figure()
#    ax1 = fig.add_subplot(211)
#    ax2 = fig.add_subplot(212)
#
#    ax1.plot(range(len(energy)), energy, **yellowline)
#    ax1.set_xlabel('N')
#    ax1.set_ylabel('energy [%s]'%value_in)
#
#    ax2.plot(range(len(energy)), energy_error, **yellowline)
#    ax2.set_xlabel('N')
#    ax2.set_ylabel('energy_error')
#    plt.show()

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

if __name__ == "__main__":
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

    options = parse_sysargs(sys.argv) 
    main(options)
