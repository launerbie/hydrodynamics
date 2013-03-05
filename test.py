#!/usr/bin/env python
import time
import numpy as np
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.gadget2.interface import Gadget2
from amuse.io.base import write_set_to_file

import plotter

def main(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
         t_end=0.01|units.day, n_steps=10):
    print N, Mtot, Rvir, t_end, n_steps
    start_time = time.time()

    converter=nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(bodies)
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + \
                hydro.thermal_energy
#??    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)
#??    write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
#??                      "hydro.hdf5", "hdf5")

    steptimes = []
    energy = [] | units.m**2*units.kg*units.s**-2
    timerange = np.linspace(0,t_end.value_in(units.day), n_steps) | units.day
    for t in timerange:
        timestep_start = time.time()
        hydro.evolve_model(t)
#??        hydro_to_framework.copy()                
#??        write_set_to_file(bodies.savepoint(t), \
#??                          "hydro.hdf5", "hdf5")

        Etot = hydro.kinetic_energy + hydro.potential_energy + \
                hydro.thermal_energy

        steptime = time.time()-timestep_start

        energy.append(Etot)
        steptimes.append(steptime)

#??        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(),
#??        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot

    hydro.stop()
    end_time = time.time()
    total_time = end_time - start_time
   
    print total_time, "seconds in simulate_system()"

    plotter.plot_steptimes(steptimes)
    plotter.plot_energy(energy)

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
    return result


if __name__ in ('__main__', '__plot__'):
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)


