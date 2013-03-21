#!/usr/bin/env python
import h5py
from amuse.units import units
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.units import quantities 

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
    """
    def __init__(self, keywords):
        """ Generates the attributes of HydroResults using the keywords 
        of the dictionary as attribute names."""
        self.__dict__.update(keywords)

class HDF5Ready(object):
    """ Wraps VectorQuantities such that they can be written by
    write_hdf5(). Basically needed to add information as attributes to
    the hdf5 files. It contains information such as the unit of the
    VectorQuantity and possibly more."""
    def __init__(self, vq, keyword):
        self.keyword = keyword
        self.nparray = vq.value_in(vq.unit)

        self.attributes = {}
        try:
            self.attributes['unit'] = vq.unit.to_reduced_form().__str__()
        except AttributeError:
            pass
        try:
            self.attributes['mass_fractions'] = vq.mf.__str__()
        except AttributeError:
            pass

def write_to_hdf5(filename, data):
    """ Writes the attibutes of HydroResults to an hdf5 file."""
    f = h5py.File(filename,'w')

    for keyword in data:
        hdf5ready = HDF5Ready(data[keyword], keyword)
        f[keyword] = hdf5ready.nparray

        for name in hdf5ready.attributes:
            f[keyword].attrs[name] = hdf5ready.attributes[name]
    f.close()
    del f

def read_from_hdf5(filename):
    """ Reads an hdf5 file and returns a HydroResult class. """
    f = h5py.File(filename,'r')

    data = {}
    for keyword in f.keys():
        unitstring = f[keyword].attrs['unit']
        evaluable_unit = parse_unitsstring(unitstring) 
        vq = f[keyword].value | eval(evaluable_unit)
        avq = AdaptingVectorQuantity()
        avq.extend(vq)

        try:
            mass_fractions = eval(f[keyword].attrs['mass_fractions'])
            setattr(avq, 'mf', mass_fractions)
        except KeyError, AttributeError:
            pass

        data.update({keyword:avq})

    f.close()
    del f
    results = HydroResults(data)
    return results

def parse_unitsstring(string):
    """ Parses the unitstring so it can be evaluated with eval to create
    a Quantity."""
    units_splitted = string.split(' * ')
    plus_units = ["units."+elem for elem in units_splitted]
    full_unitstring =  " * ".join(plus_units)
    return full_unitstring
    



