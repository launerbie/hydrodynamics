#!/usr/bin/env python

import numpy
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.fi.interface import Fi
from amuse.io.base import write_set_to_file
from amuse.io.base import read_set_from_file


Mtot=1|units.kg;Rvir=1|units.RSun
bodies = new_plummer_gas_model(N, convert_nbody=converter)
converter = nbody_system.nbody_to_si(Mtot, Rvir)
fi = Fi(converter)
fi.gas_particles.add_particles(bodies)

def main(options):
    bodies = read_set_from_file('defaults2.hdf5','hdf5') 
    converter = nbody_system.nbody_to_si(Mtot, Rvir)
    fi
    

def print_energies(hydro):
    Kh = hydro.kinetic_energy
    Th = hydro.thermal_energy
    Vh = hydro.potential_energy
    Kg = hydro.gas_particles.kinetic_energy()
    Tg = hydro.gas_particles.thermal_energy()
    Vg = hydro.gas_particles.potential_energy()
    K = hydro.particles.kinetic_energy()
    T = hydro.particles.thermal_energy()
    V = hydro.particles.potential_energy()
    total_hydro_energy = Kh+Th+Vh
    total_gas_particles_energy = Kg+Tg+Vg
    total_particles_energy = K+T+V
    print total_hydro_energy
    print total_gas_particles_energy
    print total_particles_energy



if __name__ == "__main__":
    main()





