#!/usr/bin/env python
import theme
import os
import sys
import h5py

from optparse import OptionParser

from amuse.units import units

from matplotlib import pyplot as plt
import matplotlib as mpl

from support import read_from_hdf5

def main(options):
    """ Loads HydroResults instance from an hdf5 file and plots the 
    results."""

    if options.simresults:
        filepath = options.simresults
        results = read_from_hdf5(filepath)
        filepath = new_path_extension(filepath, 'plots', 'png')
        plot_all(results, filepath = filepath)

    if options.energy_results:
        filepath = options.energy_results
        results = read_from_hdf5(filepath)
        filepath = new_path_extension(filepath, 'plots', 'png')
        plot_energy(results, filepath = filepath)
    
    if options.runtime_results:
        filepath = options.runtime_results
        results = read_from_hdf5(filepath)
        filepath = new_path_extension(filepath, 'plots', 'png')
        plot_steptimes(results, filepath = filepath)

    return 0 

def new_path_extension(filepath, newpath, extension):
    """ Returns a new path/filename derived from the old 
    path/filename."""
    filename = os.path.split(filepath)[1]
    filename = os.path.splitext(filename)[0]+"."+extension
    filepath = "/".join([newpath,filename])
    return filepath

def plot_all(results, filepath=None):
    """ Plots some attributes of a HydroResults instance and saves it
    to png."""
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
   
    fig = plt.figure(figsize=(12,12), dpi=300)

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot(times.value_in(times.unit), lr1.value_in(lr1.unit),\
             label='10%',  **theme.line)
    ax1.plot(times.value_in(times.unit), lr2.value_in(lr1.unit),\
             label='25%', **theme.blueline)
    ax1.plot(times.value_in(times.unit), lr3.value_in(lr1.unit),\
             label='50%', **theme.redline)
    ax1.plot(times.value_in(times.unit), lr4.value_in(lr1.unit),\
             label='75%', **theme.greenline)

    ax1.set_xlabel('Time in %s'%times.unit.__str__())
    ax1.set_ylabel('Lagrangian Radius in %s'%lr1.unit.__str__())
    ax1.legend(loc='best')


    ax2.plot(times.value_in(times.unit), L.value_in(Lx.unit),\
             label='L', **theme.line )
    ax2.plot(times.value_in(times.unit), Lx.value_in(Lx.unit),\
             label='Lx', **theme.blueline )
    ax2.plot(times.value_in(times.unit), Ly.value_in(Lx.unit),\
             label='Ly', **theme.redline )
    ax2.plot(times.value_in(times.unit), Lz.value_in(Lx.unit),\
             label='Lz', **theme.greenline)

    ax2.set_xlabel('Time in %s'%times.unit.__str__())
    ax2.set_ylabel('Angular Momentum in %s'%L.unit.__str__())
    ax2.legend(loc='best')


    ax3.plot(times.value_in(times.unit), Ekin.value_in(Ekin.unit),\
             label='Kinetic',  **theme.redline)
    ax3.plot(times.value_in(times.unit), Epot.value_in(Ekin.unit),\
             label='Potential', **theme.greenline)
    ax3.plot(times.value_in(times.unit), Etot.value_in(Ekin.unit),\
             label='Total',**theme.line)

    ax3.set_xlabel('Time in %s'%times.unit.__str__())
    ax3.set_ylabel('Energy in %s'%Etot.unit.__str__())
    ax3.legend(loc='best')

    ax4.plot(radius_initial.value_in(radius_initial.unit),\
             densities_initial.value_in(densities_initial.unit),\
             label='Initial', **theme.line )
    ax4.plot(radius_final.value_in(radius_final.unit),\
             densities_final.value_in(densities_initial.unit),\
             label='Final', **theme.redline )

    ax4.set_xlabel('Radius in %s'%radius_initial.unit.__str__())
    ax4.set_ylabel('Density in %s'%densities_initial.unit.__str__())
    ax4.legend(loc='best')

    plt.suptitle('Particles: %i  Steps: %i '%(N, n_steps ))

    if filepath:
        plt.savefig(filepath)

def plot_steptimes(results, filepath = None):
    """ Plot the number of particles versus the time for one
    step."""
    nr_of_particles = results.N.value_in(results.N.unit)
    steptimes = results.runtimes.value_in(units.s)
    dots = theme.dots
    dots['ms'] = 4

    fig = plt.figure(figsize=(12,12), dpi=300)
    ax1 = fig.add_subplot(111)
    ax1.plot(nr_of_particles, steptimes, **dots)
    ax1.set_xlabel('Number of particles')
    ax1.set_ylabel('Runtime in seconds for 1 step')
 
    if filepath:
        plt.savefig(filepath)

def plot_energy(results, filepath = None):
    """ Plot the number of particles versus the energy error."""
    nr_of_particles = results.N.value_in(units.no_unit)
    energy_errors = results.energy_errors.value_in(units.no_unit)

    dots = theme.dots
    dots['ms'] = 4

    fig = plt.figure(figsize=(12,12), dpi=300)
    ax1 = fig.add_subplot(111)

    ax1.plot(nr_of_particles, energy_errors, **dots)
    ax1.set_xlabel('Number of particles')
    ax1.set_ylabel('Energy error')

    if filepath:
        plt.savefig(filepath)

def parse_sysargs(sysargs):
    parser = OptionParser()
    parser.set_defaults(test=False, skip_table=False)

    parser.add_option("-f", "--file", action="store",\
                      dest="simresults", default=None,\
                      help="Simulation results HDF5 file to open")

    parser.add_option("-e", "--energy", action="store",\
                      dest="energy_results", default=None,\
                      help="Energy errors results hdf5 file")

    parser.add_option("-t", "--times", action="store", \
                      dest="runtime_results", default=None,\
                      help="Runtimes results hdf5 file")

    options, args = parser.parse_args(sysargs[1:])
    return options

if __name__ == "__main__":
    options = parse_sysargs(sys.argv) 
    main(options)
