#!/usr/bin/env python
import time
import numpy 

from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.fi.interface import Fi
from amuse.io.base import write_set_to_file

import plotter

def main(options):
    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'write_hdf5':options.write_hdf5}

    program_options = {'N_vs_t':options.N_vs_t,\
                       'N_vs_E':options.N_vs_E}

    if options.N_vs_t == True:
        create_N_vs_t()

    if options.N_vs_E == True:
        create_N_vs_E()

    run_hydrodynamics(**hydro_options)

def create_N_vs_t(print_it=False, N_list=None):
    if N_list:
        N = N_list
    else:
        #N = [10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000] 
        N = [10, 50, 100, 500, 1000, 5000] 
    total_runtimes = []
    for nr_particles in N: 
        start_time = time.time()
        run_hydrodynamics(N=nr_particles, n_steps=1)
        total_runtime = time.time() - start_time
        total_runtimes.append((total_runtime, nr_particles))
    if print_it == True:
        for elem in total_runtimes:
            print "runtime:", elem[0],"nr_particles", elem[1]

def create_N_vs_E(print_it=False, N_list=None):
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

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
         t_end=0.5|units.day, n_steps=10, write_hdf5=False):
    print "N:%i, Mtot:%s, Rvir:%s, t_end:%s, n_steps:%s"%(N, \
           Mtot.as_string_in(units.MSun), \
           Rvir.as_string_in(units.RSun), \
           t_end.as_string_in(units.day), n_steps)

    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    hydro = Fi(converter)
    hydro.gas_particles.add_particles(bodies)
#    hydro.parameters.isothermal_flag = True
#    hydro.parameters.integrate_entropy_flag = False
#    hydro.parameters.gamma = 1

    if write_hdf5 == True:
       hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)
       write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
                         "hydro.hdf5", "hdf5")

       timerange = numpy.linspace(0, t_end.value_in(units.day),\
                                  n_steps) | units.day
       for t in timerange:
           hydro.evolve_model(t)
           hydro_to_framework.copy()                
           write_set_to_file(bodies.savepoint(t),"hydro.hdf5", "hdf5")
       hydro.stop()

    #energy_error = 1.0-(Etot_end/Etot_init)
    return 0 

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int", default = 1000,
                      help="number of stars [1000]")
    result.add_option("-n", dest="n_steps", type="int", default = 100,
                      help="number of steps [200]")
    result.add_option("-t", unit=units.day,
                      dest="t_end", type="float", default = 0.5|units.day,
                      help="end time of the simulation [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mtot", type="float", default = 1|units.MSun,
                      help="Mass of molcular cloud [%default]")
    result.add_option("-R", unit=units.RSun,
                      dest="Rvir", type="float", default = 1|units.RSun,
                      help="Radius of cloud [%default]")
    result.add_option("-T", dest="N_vs_t", action='store_true',\
                      default = False,\
                      help="Creates plot of one-step time vs N.")
    result.add_option("-E", dest="N_vs_E", action='store_true',\
                      default = False,\
                      help="Creates plot of N vs energy error.")
    result.add_option("-H", dest="write_hdf5", action='store_true',\
                      default = False,\
                      help="Writes results to an hdf5 file.")
    return result

if __name__ in ('__main__', '__plot__'):
    options, arguments = new_option_parser().parse_args()
    main(options)


#    Etot_init = hydro.kinetic_energy + hydro.potential_energy + \
#                hydro.thermal_energy

#    Etot_end = hydro.kinetic_energy + hydro.potential_energy + \
#            hydro.thermal_energy
#        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(),
#        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot


