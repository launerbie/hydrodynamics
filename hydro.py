#!/usr/bin/env python
import numpy
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.fi.interface import Fi
from amuse.io.base import write_set_to_file

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                      t_end=0.5|units.day, n_steps=10, write_hdf5=None):

    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    fi = Fi(converter)
    fi.gas_particles.add_particles(bodies)
#    fi.parameters.isothermal_flag = True
#    fi.parameters.integrate_entropy_flag = False
#    fi.parameters.gamma = 1

    timerange = numpy.linspace(0, t_end.value_in(t_end.unit),\
                                  n_steps) | t_end.unit

    if write_hdf5:
       filename = write_hdf5
       fi_to_framework = fi.gas_particles.new_channel_to(bodies)
       write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
                         filename, "hdf5")
       for t in timerange:
           #print "Evolving to t=%s"%t.as_string_in(t.unit)
           fi.evolve_model(t)
           fi_to_framework.copy()
           write_set_to_file(bodies.savepoint(t), filename, "hdf5")
    else:
       for t in timerange:
           #print "Evolving to t=%s"%t.as_string_in(t.unit)
           fi.evolve_model(t)

    fi.stop()

    #energy_error = 1.0-(Etot_end/Etot_init)
    return 0

