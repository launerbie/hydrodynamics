#!/usr/bin/env python
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

from samples import init_conditions

#q=VectorQuantity.new_from_scalar_quantities(*b)

class HydroResults(object):
    def __init__(self, results):
        self.angular_momentum = results['angular_momentum']
        self.potential_energy= results['potential_energy']
        self.kinetic_energy = results['kinetic_energy'] 
        self.total_energy = results['total_energy']
        self.radial_profile = "think of something"
        self.lagrangianradii = results['lagrangianradii']

    def write_to_hdf5(self, filename):
        """ Write all of the above to hdf5 file."""
        f = h5py.File(filename,'w')

        f['angular_momentum'] = self.angular_momentum.value_in(\
                                    self.angular_momentum.unit)
        f['angular_momentum'].attrs['unit'] = 'someunit'

        f['potential_energy'] = self.potential_energy.value_in(\
                                    self.potential_energy.unit)
        f['potential_energy'].attrs['unit'] = 'someunit_energy'

        f['kinetic_energy'] = self.kinetic_energy.value_in(\
                                    self.kinetic_energy.unit)
        f['kinetic_energy'].attrs['unit'] = 'someunit_energy'

        f['total_energy'] = self.total_energy.value_in(\
                                    self.total_energy.unit)
        f['total_energy'].attrs['unit'] = 'someunit_energy'

        for i, radius in enumerate(self.lagrangianradii[0][1]):
            key ='lagrangianradii'+"_"+str(radius)
            print "key:",key 
            f[key] = self.lagrangianradii
            f[key].attrs['unit'] = 'meter'
        f.close()
        del f

    def read_from_hdf5(self, filename):
        """ Read all of the above from file. """
        pass

class LagrangianRadii(object):
      def __init__(self, list_):
          radii = list_[0][0]
          lr = [elem[0] for elem in list_]
          for i, radius in radii:
              pass
          pass

def main(options):
    start_time = time.time()
    hydro_options = {'N':options.N, 'Mtot':options.Mtot,\
                     'Rvir':options.Rvir, 't_end':options.t_end,\
                     'n_steps':options.n_steps,\
                     'write_hdf5':options.write_hdf5}

    if options.N_vs_t == True:
        create_N_vs_t()

    if options.N_vs_E == True:
        create_N_vs_E()

    data, r = run_hydrodynamics(**hydro_options)
    
    end_time = time.time()
    print "Total runtime:", end_time-start_time, "seconds."
    
    r.write_to_hdf5(options.results_out)

    return 0

def run_hydrodynamics(N=100, Mtot=1|units.MSun, Rvir=1|units.RSun,
                      t_end=0.5|units.day, n_steps=10, write_hdf5=None):

    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    bodies = new_plummer_gas_model(N, convert_nbody=converter)

    fi = Fi(converter)
    fi.gas_particles.add_particles(bodies)

    #Adiabetic equation of state means the following?:
    fi.parameters.isothermal_flag = True
    fi.parameters.integrate_entropy_flag = False
    #fi.parameters.gamma = 1

    data = {'lagrangianradii':[],\
            'angular_momentum':AdaptingVectorQuantity(),\
            'radial_profile_initial':AdaptingVectorQuantity(),\
            'radial_profile_final':AdaptingVectorQuantity(),\
            'kinetic_energy':AdaptingVectorQuantity(),\
            'potential_energy':AdaptingVectorQuantity(),\
            'total_energy':AdaptingVectorQuantity() } 

    data['radial_profile_initial'].append(1|units.m)
    data['radial_profile_final'].append(1|units.m)

    mass_fraction = [0.10, 0.25, 0.50, 0.75]

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
                                          mf=mass_fraction))
           
    fi.stop()
    #energy_error = 1.0-(Etot_end/Etot_init)
    results=data
    return results, HydroResults(data)

def create_N_vs_t(print_it=False, N_list=None):
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
    #writeto_hdf5file

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
    #writeto_hdf5file

def new_option_parser():
    from amuse.units.optparse import OptionParser
    from amuse.units import units
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
    result.add_option("-H", dest="write_hdf5", action='store',\
                      default = None,\
                      help="Specifies filename for hdf5 output.")
    result.add_option("-P", dest="results_out", action='store',\
                      default = None,\
                      help="Specifies filename for plotting data (hdf5 format).")
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


#           print type(fi.gas_particles.total_angular_momentum())
#           print type(fi.kinetic_energy)
#           print type(fi.potential_energy)
#           print type(fi.total_energy)
#           print type(fi.particles.LagrangianRadii(unit_converter=converter,\
#                     mf=[0.10, 0.25, 0.50, 0.75]))
