#!/usr/bin/env python

import os
import time
import numpy
from progressbar import progressbar as pb

from support import HydroResults, write_to_hdf5
from support import radial_density
from support import setup_directories
from viewer import plot_all
from viewer import plot_steptimes
from viewer import plot_energy

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.io.base import write_set_to_file, read_set_from_file
from amuse.community.fi.interface import Fi

def main(options):
    """ Separates the arguments for the hydro_solver from the
    program control-flow arguments. Runs stuff."""

    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'bodyname':options.bodyname,\
                     'vx':options.vx,\
                     'vy':options.vy,\
                     'vz':options.vz,\
                     'plummer1':options.plummer1,\
                     'plummer2':options.plummer2}

    setup_directories('bodies','hydroresults','plots')

    if options.N_vs_t:
        create_N_vs_t(options.N_vs_t)

    elif options.N_vs_E:
        create_N_vs_E(options.N_vs_E)

    else:
        results = run_hydrodynamics(**hydro_options)

    if options.results_out:
        if os.path.dirname(options.results_out) == '':
            filename ='hydroresults/'+options.results_out
        write_to_hdf5(filename, results.__dict__)
        pngfilename = options.results_out.split('.')[0]+'.png'
        plot_all(results, filepath = "plots/"+pngfilename)
    return 0

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                      t_end=0.5|units.day, n_steps=10,\
                      vx = 0 |(units.RSun/units.day),\
                      vy = 0 |(units.RSun/units.day),\
                      vz = 0 |(units.RSun/units.day),\
                      plummer1=None, plummer2=None,\
                      bodyname = None):

    """ Runs the hydrodynamics simulation and returns a HydroResults
    instance. 


    FUNCTION WALKTHROUGH:

    In the following explanation 'plummer1' and 'plummer2' are assumed
    to be hdf5 files written by the function write_set_to_file().

    Case 1: 
    If 'plummer1' and 'plummer2' are filenames of hdf5 files, then these 
    two plummer spheres will be smashed together. 

    Case 2:
    If only plummer1 is supplied, then it will evolve plummer1 with 
    t_end timerange and n_steps steps.

    Case 3:
    If no plummers spheres are supplied, then it will generate a new
    plummer sphere using the default/supplied initial conditions. 

   
    OUTPUT FILES:

    If 'results_out' is specified, the HydroResult instance is written
    to file in HDF5 format. This however does not use 
    write_set_to_file() which writes the entire Particles class and
    its attributes to file at each dt, but uses write_to_hdf5() 
    from the 'support' module which is tailored to writing
    HydroResults instances to file. This HDF5 contains all necessary
    data to plot the required plots of the assignment.
    
    In addition, the last snapshot of the Particles instance is written
    to file using write_set_to_file(), the latter file is written to 
    the 'bodies' directory. Only the last snapshot is written to file
    because the history of a Particle set is not of interest when 
    reloading them to smash plummer spheres.   """

    converter = nbody_system.nbody_to_si(Mtot, Rvir)

    fi = Fi(converter)

    if plummer1 and plummer2: 
        eta_smash = 0.3 |units.day
        if plummer1 == plummer2:
            bodies1 = read_set_from_file(plummer1, format='hdf5')
            bodies2 = bodies1.copy()
            bodies2.key += 1
        else:
            bodies1 = read_set_from_file(plummer1, format='hdf5')
            bodies2 = read_set_from_file(plummer2, format='hdf5')

        bodies1.move_to_center()
        bodies2.move_to_center()

        if vx.value_in(vx.unit) == 0  and vy.value_in(vy.unit) == 0 \
            and vz.value_in(vz.unit) == 0:
            bodies1.x += -1 |units.RSun
            bodies2.x += 1 |units.RSun

        else:
            bodies1.x += (-1)*vx*eta_smash
            bodies2.x += 1*vx*eta_smash
            bodies1.vx += vx
            bodies2.vx += (-1)*vx 
            bodies1.vy += vy
            bodies2.vy += (-1)*vy
            bodies1.vz += vz
            bodies2.vz += (-1)*vz

        bodies1.add_particles(bodies2)
        bodies = bodies1

    elif plummer1 or plummer2: 
        if plummer1:
            bodies = read_set_from_file(plummer1)
        else:
            bodies = read_set_from_file(plummer2)
        bodies.move_to_center()

    else: 
        bodies = new_plummer_gas_model(N, convert_nbody=converter)

    fi.gas_particles.add_particles(bodies)
    fi_to_framework = fi.gas_particles.new_channel_to(bodies)

    fi.parameters.self_gravity_flag = True

    data = {'lagrangianradii':AdaptingVectorQuantity(),\
            'angular_momentum':AdaptingVectorQuantity(),\
            'time':AdaptingVectorQuantity(),\
            'positions':AdaptingVectorQuantity(),\
            'kinetic_energy':AdaptingVectorQuantity(),\
            'potential_energy':AdaptingVectorQuantity(),\
            'total_energy':AdaptingVectorQuantity() } 

    mass_fractions = [0.10, 0.25, 0.50, 0.75]
    setattr(data['lagrangianradii'], 'mf', mass_fractions)

    data['radius_initial'], data['densities_initial'] = radial_density(\
         fi.particles.x, fi.particles.mass, N=10, dim=1)

    timerange = numpy.linspace(0, t_end.value_in(t_end.unit),\
                                  n_steps) | t_end.unit
    data['time'].extend(timerange)
 
    fi.parameters.timestep = t_end/(n_steps+1)

    widget = drawwidget("Evolving")
    pbar = pb.ProgressBar(widgets=widget, maxval=len(timerange)).start()

    for i, t in enumerate(timerange):
        fi.evolve_model(t)
        data['kinetic_energy'].append(fi.kinetic_energy)
        data['potential_energy'].append(fi.potential_energy)
        data['total_energy'].append(fi.total_energy)
        data['positions'].append(fi.particles.position)
        data['angular_momentum'].append(fi.gas_particles.\
                                        total_angular_momentum())
        data['lagrangianradii'].append(fi.particles.LagrangianRadii(\
                                       unit_converter=converter,\
                                       mf=mass_fractions)[0])
        if t == timerange[-1] and bodyname:
            if os.path.dirname(bodyname) == '':
                filename = "bodies/"+bodyname 
            fi_to_framework.copy()
            write_set_to_file(bodies.savepoint(t), filename, "hdf5")

        pbar.update(i)
    pbar.finish()

    data['radius_final'], data['densities_final'] = radial_density(\
         fi.particles.x, fi.particles.mass, N=10, dim=1)

    fi.stop()

    results = HydroResults(data)
    return results 

def create_N_vs_t(filename, N_list=None):
    """ Iteratively runs run_hydrodynamics for several N in order to 
    create the data for the N vs. wallclock-time plot. Writes data to 
    hdf5 file."""

    if N_list:
        N = N_list
    else:
        a = range(100,1000,100)
        b = range(1000,10000,1000)
        N = a+b

    total_runtimes = AdaptingVectorQuantity() 
    for nr_particles in N: 
        start_time = time.time() |units.s
        run_hydrodynamics(N=nr_particles, n_steps=1)
        end_time = time.time() |units.s
        total_runtime = end_time - start_time
        total_runtimes.append(total_runtime)

    N_as_vq = N | units.no_unit 

    data = {'N':N_as_vq, 'runtimes':total_runtimes}
    results = HydroResults(data)

    write_to_hdf5('hydroresults/'+filename,\
                  results.__dict__)
    pngfilename = options.N_vs_t.split('.')[0]+'.png'
    plot_steptimes(results, filepath = "plots/"+pngfilename)

def create_N_vs_E(filename, N_list = None):
    """ Iteratively runs run_hydrodynamics for several N in order to 
    create the data for the N vs. Energy error plot. Writes data to hdf5
    file."""

    if N_list:
        N = N_list
    else:
        N = range(100, 1000, 50) 

    energy_errors = AdaptingVectorQuantity() 
    for nr_particles in N: 
        results = run_hydrodynamics(N=nr_particles, n_steps=50)
        energy_begin = results.total_energy[0]
        energy_end = results.total_energy[-1]
        energy_error = (energy_begin-energy_end)/energy_end
        energy_error = energy_error |units.no_unit
        energy_errors.append(energy_error)

    N_as_vq = N | units.no_unit 

    data = {'N':N_as_vq, 'energy_errors':energy_errors}
    results = HydroResults(data)

    write_to_hdf5('hydroresults/'+filename,\
                  results.__dict__)
    pngfilename = options.N_vs_E.split('.')[0]+'.png'
    plot_energy(results, filepath = "plots/"+pngfilename)

def drawwidget(proces_discription):
    widgets = [proces_discription+": ", pb.Percentage(), ' ',
               pb.Bar(marker='#',left='[',right=']'),
               ' ', pb.ETA()]
    return widgets

def new_option_parser():
    from amuse.units.optparse import OptionParser
    from amuse.units import units
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 1000,
                      help="number of stars [1000]")
    result.add_option("-n", dest="n_steps", type="int", default = 100,
                      help="number of steps [200]")
    result.add_option("-t", unit=units.day,
                      dest="t_end", type="float", \
                      default = 0.5|units.day,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of molcular cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 1|units.RSun,
                      help="Radius of cloud [%default]")
    result.add_option("-T", dest="N_vs_t", action='store',\
                      default = None, type='string',\
                      help="Filename for N vs time hdf5 file.")
    result.add_option("-E", dest="N_vs_E", action='store',\
                      default = None, type='string',\
                      help="Filename for N vs Energy hdf5 file.")
    result.add_option("-H", dest="results_out", action='store',\
                      default = None,\
                      help="path/to/filename for plotting data (hdf5 format).")
    result.add_option("-B", dest="bodyname", action='store',\
                      default = None,\
                      help="Filename for resulting body (hdf5 format).")
    result.add_option("-p", dest="plummer1",type="string", \
                      action='store', default = None,\
                      help="First plummer sphere.")
    result.add_option("-q", dest="plummer2",type="string", \
                      action='store', default = None,\
                      help="Second plummer sphere.")
    result.add_option("--vx", unit=(units.RSun/units.day),
                      dest="vx", type="float", \
                      default = 0|(units.RSun/units.day),
                      help="Smash velocity x in RSun/day.[%default]")
    result.add_option("--vy", unit=(units.RSun/units.day),
                      dest="vy", type="float", \
                      default = 0|(units.RSun/units.day),
                      help="Smash velocity y in RSun/day.[%default]")
    result.add_option("--vz", unit=(units.RSun/units.day),
                      dest="vz", type="float", \
                      default = 0|(units.RSun/units.day),
                      help="Smash velocity z in RSun/day.[%default]")
    return result

if __name__ in ('__main__'):
    options, arguments = new_option_parser().parse_args()
    main(options)
