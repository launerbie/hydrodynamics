#!/usr/bin/env python
#This line is 72 characters loooooooooooooooooooooooooooooooooooooooong.
import numpy
import time
from support import HydroResults, write_to_hdf5

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units.quantities import VectorQuantity
from amuse.io.base import write_set_to_file
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.radial_profile import radial_density
from amuse.community.fi.interface import Fi

def main(options):
    """ Separates the arguments for the hydro_solver from the
    program control-flow arguments. Runs stuff."""
    start_time = time.time()
    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'write_hdf5':options.write_hdf5}

    results = smash_plummers(**hydro_options)

    end_time = time.time()
    print "Total runtime:", end_time-start_time, "seconds."

    if options.results_out:
        write_to_hdf5(options.results_out, results.__dict__)
    return 0

def smash_plummers(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                   t_end=0.5|units.day, n_steps=10, write_hdf5=None,\
                   plummer1=None, plummer2=None):
    """Runs the hydrodynamics simulation and returns a HydroResults
    instance."""

    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    fi = Fi(converter)

    #if plummer1 and no plummer2
    # run as usual with just plummer1

    #if no plummer1 and no plummer2
    # generate 1 plummer like before from initial conditions
    # passed to function

    #if no plummer1 and yes plummer2
    # run as usual with plummer2 as primary body

    #if plummer1 and plummer2
    # Smash them.

    bodies = new_plummer_gas_model(N, convert_nbody=converter)
    bodies2 = new_plummer_gas_model(N, convert_nbody=converter)

    bodies2.x += 2 | units.RSun
    bodies2.vx = -0.0 | (units.RSun / units.day)

    fi.gas_particles.add_particles(bodies)
    fi.gas_particles.add_particles(bodies2)
   
    print fi.gas_particles

    fi.parameters.self_gravity_flag = True

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

    if write_hdf5:
       filename = write_hdf5
       fi_to_framework = fi.gas_particles.new_channel_to(bodies)
       write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
                         filename, "hdf5")
       for t in timerange:
           print "Evolving to t=%s"%t.as_string_in(t.unit)
           fi.evolve_model(t)
           fi_to_framework.copy()
           write_set_to_file(bodies.savepoint(t), filename, "hdf5")
    else:
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
    result.add_option("-H", dest="write_hdf5", type="string", 
                      action='store', default = None,\
                      help="Filename for hdf5 output.")
    result.add_option("-P", dest="results_out", action='store',\
                      default = None,\
                      help="Filename for plotting data (hdf5 format).")
#    result.add_option("-p1", dest="plummer1", action='store',\
#                      default = None,\
#                      help="HDF5 file for plummer sphere one.")
#    result.add_option("-p2", dest="plummer2", action='store',\
#                      default = None,\
#                      help="HDF5 file for plummer sphere two.")
    return result

if __name__ in ('__main__'):
    options, arguments = new_option_parser().parse_args()
    main(options)

