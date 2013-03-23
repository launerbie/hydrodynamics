#!/usr/bin/env python
#This line is 72 characters loooooooooooooooooooooooooooooooooooooooong.

import time
import numpy
from support import HydroResults, write_to_hdf5

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.ext.radial_profile import radial_density
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.io.base import write_set_to_file, read_set_from_file
from amuse.community.fi.interface import Fi
#todo create runset by datetimestamp

def main(options):
    """ Separates the arguments for the hydro_solver from the
    program control-flow arguments. Runs stuff."""
    start_time = time.time()
    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'write_hdf5':options.write_hdf5,\
                     'smashvelocity':options.smashvelocity,\
                     'plummer1':options.plummer1,\
                     'plummer2':options.plummer2}

    if options.N_vs_t:
        create_N_vs_t(options.N_vs_t)

    if options.N_vs_E:
        create_N_vs_E(options.N_vs_E)

    results = run_hydrodynamics(**hydro_options)
    
    end_time = time.time()
    print "Total runtime:", end_time-start_time, "seconds."
    
    if options.results_out:
        write_to_hdf5(options.results_out, results.__dict__)
    return 0

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                      t_end=0.5|units.day, n_steps=10,\
                      smashvelocity = 2 |(units.RSun/units.day),\
                      write_hdf5=None, plummer1=None, plummer2=None):
    """ Runs the hydrodynamics simulation and returns a HydroResults
    instance. 

    Function walkthrough:

    In this function 'plummer1' and 'plummer2' are assumed to be hdf5 
    files written by the function write_set_to_file().

    Case 1: 
    If 'plummer1' and 'plummer2' are filenames of hdf5 files, then these 
    two plummer spheres will be smashed together. 

    Case 2:
    If only plummer1 is supplied, then it will evolve plummer1 with 
    t_end timerange and n_steps steps.

    Case 3:
    If no plummers spheres are supplied, then it will generate a new
    plummer sphere using the default/supplied initial conditions. 

   
    OUTPUT files:

    Case 1:
    If write_hdf5 is a specified output filename, then it will create
    an hdf5 file which contains snapshots of the particle set at each
    timestep.

    Case 2:
    If 
     
    """

    converter = nbody_system.nbody_to_si(Mtot, Rvir)

    fi = Fi(converter)

    if plummer1 and plummer2: 
        eta_smash = 0.3 |units.day
        bodies1 = read_set_from_file(plummer1, format='hdf5')
        bodies2 = read_set_from_file(plummer2, format='hdf5')
        bodies1.move_to_center()
        bodies2.move_to_center()
        bodies1.x += (-1)*smashvelocity*eta_smash
        bodies2.x += 1*smashvelocity*eta_smash
        bodies1.vx += smashvelocity
        bodies2.vx += (-1)*smashvelocity 
        bodies1.vy += smashvelocity*0.3
        bodies2.vy += (-1)*smashvelocity *0.3
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
    #fi.parameters.stopping_condition_minimum_density = 0.0|(units.kg*units.m**-3)

    data = {'lagrangianradii':AdaptingVectorQuantity(),\
            'angular_momentum':AdaptingVectorQuantity(),\
            'time':AdaptingVectorQuantity(),\
            'positions':AdaptingVectorQuantity(),\
            'kinetic_energy':AdaptingVectorQuantity(),\
            'potential_energy':AdaptingVectorQuantity(),\
            'total_energy':AdaptingVectorQuantity() } 

    mass_fractions = [0.10, 0.25, 0.50, 0.75]

    radius_init, densities_init = radial_density(fi.particles.x, \
                                      fi.particles.mass, N=10, dim=1)

    radius_init = VectorQuantity.new_from_scalar_quantities(*radius_init)
    densities_init = VectorQuantity.new_from_scalar_quantities(*densities_init)
    
    data['radius_initial'] = radius_init 
    data['densities_initial'] = densities_init 

    timerange = numpy.linspace(0, t_end.value_in(t_end.unit),\
                                  n_steps) | t_end.unit
    data['time'].extend(timerange)
 
    fi.parameters.timestep = t_end/(n_steps+1)
    
    if write_hdf5:
       filename = write_hdf5
       write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
                         filename, "hdf5")
       for t in timerange:
           print "Evolving to t=%s"%t.as_string_in(t.unit)
           fi.evolve_model(t)
           fi_to_framework.copy()
           write_set_to_file(bodies.savepoint(t), filename, "hdf5")
    else:
       filename = "body_N"+str(N)+"n"+str(n_steps)+".hdf5" 
       for t in timerange:
           print "Evolving to t=%s"%t.as_string_in(t.unit)
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
           if t == timerange[-1]:
               fi_to_framework.copy()
               write_set_to_file(bodies.savepoint(t), filename, "hdf5")
           
    setattr(data['lagrangianradii'], 'mf', mass_fractions)

    radius_end, densities_end = radial_density(fi.particles.x, \
                                      fi.particles.mass, N=10, dim=1)

    radius_end = VectorQuantity.new_from_scalar_quantities(*radius_end)
    densities_end = VectorQuantity.new_from_scalar_quantities(*densities_end)
    
    data['radius_final'] = radius_end 
    data['densities_final'] = densities_end 

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
        N = [10, 50, 100, 500, 1000, 5000, 10000] 
    total_runtimes = []
    for nr_particles in N: 
        start_time = time.time()
        run_hydrodynamics(N=nr_particles, n_steps=1)
        total_runtime = time.time() - start_time
        total_runtimes.append((total_runtime, nr_particles))
    if print_it == True:
        for elem in total_runtimes:
            print "runtime:", elem[0],"nr_particles", elem[1]

    N_as_vq = N | units.no_unit 
    total_runtimes_as_vq = total_runtimes | units.s

    data = {'N':N_as_vq, 'runtimes':total_runtimes_as_vq}
    results = HydroResults(data)
    results.writeto_hdf5file(filename)

def create_N_vs_E(filename, N_list=None):
    """ Iteratively runs run_hydrodynamics for several N in order to 
    create the data for the N vs. Energy error plot. Writes data to hdf5
    file."""
    if N_list:
        N = N_list
    else:
        N = range(10, 1000, 50) 
    energy_errors = []
    for nr_particles in N: 
        energy_error = run_hydrodynamics(N=nr_particles, n_steps=100)
        energy_errors.append((energy_error, nr_particles))

    if print_it == True:
        for elem in energy_errors:
            print "energy error:", elem[0],"nr_particles", elem[1]

    N_as_vq = N | units.no_unit 
    energy_errors_as_vq = energy_errors | units.s

    data = {'N':N_as_vq, 'energyerror':energy_errors_as_vq}
    results = HydroResults(data)
    results.writeto_hdf5file(filename)

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
                      default = None,\
                      help="Create plot of one-step time vs N.")
    result.add_option("-E", dest="N_vs_E", action='store',\
                      default = None,\
                      help="Create plot of N vs energy error.")
    result.add_option("-H", dest="write_hdf5", type="string", action='store',\
                      default = None,\
                      help="Filename for hdf5 output.")
    result.add_option("-P", dest="results_out", action='store',\
                      default = None,\
                      help="Filename for plotting data (hdf5 format).")
    result.add_option("-p", dest="plummer1",type="string", action='store',\
                      default = None,\
                      help="First plummer sphere.")
    result.add_option("-q", dest="plummer2",type="string", action='store',\
                      default = None,\
                      help="Second plummer sphere.")
    result.add_option("-v", unit=(units.RSun/units.day),
                      dest="smashvelocity", type="float", \
                      default = 2|(units.RSun/units.day),
                      help="Smash velocity in RSun/day.[%default]")
    return result

if __name__ in ('__main__'):
    options, arguments = new_option_parser().parse_args()
    main(options)

    #fi.parameters.isothermal_flag = True
    #fi.parameters.integrate_entropy_flag = False
    #fi.parameters.gamma = 1
