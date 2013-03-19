#!/usr/bin/env python
#This line is 72 characters loooooooooooooooooooooooooooooooooooooooong.
import h5py
import time
import numpy
import plotter

from amuse.units import units
from amuse.units import nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units.quantities import VectorQuantity
from amuse.io.base import write_set_to_file
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.radial_profile import radial_density
from amuse.community.fi.interface import Fi

class HydroResults(object):
    """ Container for the results of the hydrodynamics simulation. 
    Instantiate a HydroResult class by passing it a dictionary of
    VectorQuantities.
 
    Example usage:
    >>> from amuse.lab import *
    >>> from hydro_fi import HydroResults
    >>> from hydro_fi import HDF5Ready

    >>> data = dict(masses=[1,2,3]|units.kg, positions=[1,2,3]|units.m)
    >>> results = HydroResults(data)
    >>> results.masses
    quantity<[1, 2, 3] kg>
    >>> results.positions
    quantity<[1, 2, 3] m>
    >>> results.write_to_hdf5('testdata.hdf5')

    """

    def __init__(self, keywords):
        """ Generates the attributes of HydroResults using the keywords 
        of the dictionary as attribute names."""
        self.__dict__.update(keywords)

    def write_to_hdf5(self, filename):
        """ Writes the attibutes of HydroResults to an hdf5 file."""
        f = h5py.File(filename,'w')

        for keyword in self.__dict__:
            hdf5ready = HDF5Ready(self.__dict__[keyword], keyword)
            f[keyword] = hdf5ready.nparray 

            for name in hdf5ready.attributes:
                f[keyword].attrs[name] = hdf5ready.attributes[name]
        f.close()
        del f

class HDF5Ready(object):
    """ Wraps VectorQuantities such that they can be written by
    write_hdf5(). Basically needed to add information as attributes to
    the hdf5 files. It contains information such as the unit of the
    VectorQuantity."""
    def __init__(self, vq, keyword):
        self.keyword = keyword
        self.nparray = vq.value_in(vq.unit) 

        self.attributes = {} 
        try:
            self.attributes['unit'] = vq.unit.__str__() 
        except AttributeError:
            pass
        try:
            self.attributes['mass_fractions'] = vq.mf.__str__() 
        except AttributeError:
            pass

def main(options):
    """ Separates the arguments for the hydro_solver from the
    program control-flow arguments. Runs stuff."""
    start_time = time.time()
    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'write_hdf5':options.write_hdf5}

    if options.N_vs_t == True:
        create_N_vs_t()

    if options.N_vs_E == True:
        create_N_vs_E()

    results = run_hydrodynamics(**hydro_options)
    
    end_time = time.time()
    print "Total runtime:", end_time-start_time, "seconds."
    
    results.write_to_hdf5(options.results_out)

    return 0

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                      t_end=0.5|units.day, n_steps=10, write_hdf5=None):
    """Runs the hydrodynamics simulation and returns a HydroResults
    instance."""

    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    fi = Fi(converter)
    fi.gas_particles.add_particles(bodies)

    #Adiabetic equation of state means the following?:
    fi.parameters.isothermal_flag = True
    fi.parameters.integrate_entropy_flag = False
    #fi.parameters.gamma = 1

    data = {'lagrangianradii':AdaptingVectorQuantity(),\
            'angular_momentum':AdaptingVectorQuantity(),\
            'radial_profile_initial':AdaptingVectorQuantity(),\
            'radial_profile_final':AdaptingVectorQuantity(),\
            'kinetic_energy':AdaptingVectorQuantity(),\
            'potential_energy':AdaptingVectorQuantity(),\
            'total_energy':AdaptingVectorQuantity() } 

    data['radial_profile_initial'].append(1|units.m) # placeholder
    data['radial_profile_final'].append(1|units.m) # placeholder
    mass_fractions = [0.10, 0.25, 0.50, 0.75]

    timerange = numpy.linspace(0, t_end.value_in(t_end.unit),\
                                  n_steps) | t_end.unit

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
           data['angular_momentum'].append(fi.gas_particles.\
                                           total_angular_momentum())
           data['lagrangianradii'].append(fi.particles.LagrangianRadii(\
                                          unit_converter=converter,\
                                          mf=mass_fractions)[0])
           
    fi.stop()
    #energy_error = 1.0-(Etot_end/Etot_init)
    setattr(data['lagrangianradii'], 'mf', mass_fractions)
    results = HydroResults(data)
    return results 

def create_N_vs_t(print_it=False, N_list=None):
    """ Iteratively runs run_hydrodynamics for several N in order to 
    create the data for the N vs. wallclock-time plot. Writes data to 
    hdf5 file."""
    if N_list:
        N = N_list
    else:
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
    #writeto_hdf5file? or just plot and leave?

def create_N_vs_E(print_it=False, N_list=None):
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
    #writeto_hdf5file? or just plot and leave?

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
    result.add_option("-T", dest="N_vs_t", action='store_true',\
                      default = False,\
                      help="Create plot of one-step time vs N.")
    result.add_option("-E", dest="N_vs_E", action='store_true',\
                      default = False,\
                      help="Create plot of N vs energy error.")
    result.add_option("-H", dest="write_hdf5", action='store',\
                      default = None,\
                      help="Filename for hdf5 output.")
    result.add_option("-P", dest="results_out", action='store',\
                      default = None,\
                      help="Filename for plotting data (hdf5 format).")
    return result

if __name__ in ('__main__'):
    options, arguments = new_option_parser().parse_args()
    main(options)

#    Etot_init = hydro.kinetic_energy + hydro.potential_energy + \
#                hydro.thermal_energy

#    Etot_end = hydro.kinetic_energy + hydro.potential_energy + \
#            hydro.thermal_energy
#        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(),
#        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot


#from example import conditions
#q=VectorQuantity.new_from_scalar_quantities(*b)
