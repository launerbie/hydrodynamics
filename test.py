#!/usr/bin/env python
import time
import numpy as np

from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.fi.interface import Fi
from amuse.io.base import write_set_to_file

import plotter

def main(options):
    if options.__dict__['one_step'] == True:

        #N = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000] 
        N = [10, 50, 100, 500, 1000] 
        total_runtimes = []
        for nr_particles in N: 
            start_time = time.time()
            run_hydrodynamics(N=nr_particles, n_steps=1)
            total_runtime = time.time() - start_time
            total_runtimes.append((total_runtime, nr_particles))

        for elem in total_runtimes:
            print "runtime:", elem[0],"nr_particles", elem[1]

    if options.__dict__['energy_error'] == True:
        N = range(10,1000,10) 
        energy_errors = []
        for nr_particles in N: 
            energy_error = run_hydrodynamics(N=nr_particles, n_steps=100)
            energy_errors.append((energy_error, nr_particles))

        for elem in energy_errors:
            print "energy error:", elem[0],"nr_particles", elem[1]
        
    #plotter(total_runtimes)
    #plotter(energy_errors)

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
         t_end=0.01|units.day, n_steps=10):
    print "N:%i, Mtot:%s, Rvir:%s, t_end:%s, n_steps:%s"%(N, \
           Mtot.as_string_in(units.MSun), \
           Rvir.as_string_in(units.RSun),\
           t_end.as_string_in(units.day), n_steps)

    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    hydro = Fi(converter)
    hydro.gas_particles.add_particles(bodies)
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.isothermal_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.gamma = 1

    Etot_init = hydro.kinetic_energy + hydro.potential_energy + \
                hydro.thermal_energy

#   hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)
#   write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
#                     "hydro.hdf5", "hdf5")

    timerange = np.linspace(0,t_end.value_in(units.day), n_steps) | units.day
    for t in timerange:
        hydro.evolve_model(t)
#        hydro_to_framework.copy()                
#        write_set_to_file(bodies.savepoint(t), \
#                          "hydro.hdf5", "hdf5")

    Etot_end = hydro.kinetic_energy + hydro.potential_energy + \
            hydro.thermal_energy
#        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(),
#        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot
    hydro.stop()

    energy_error = 1.0-(Etot_end/Etot_init)
    return energy_error 

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 1000,
                      help="number of stars [1000]")
    result.add_option("-n", dest="n_steps", type="int", default = 200,
                      help="number of steps [200]")
    result.add_option("-t", unit=units.day,
                      dest="t_end", type="float", default = 0.01|units.day,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of molcular cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 1|units.RSun,
                      help="Radius of cloud [%default]")
    result.add_option("-O", dest="one_step", action='store_true',\
                      default = False,\
                      help="Just runs one integration step.")
    result.add_option("-E", dest="energy_error", action='store_true',\
                      default = False,\
                      help="Just runs one integration step.")
    return result

if __name__ in ('__main__', '__plot__'):
    options, arguments = new_option_parser().parse_args()
    main(options)


