#!/usr/bin/env python
import os
import sys
import subprocess
from optparse import OptionParser

def main(options):
    """ Generates some example results and animations.
    Must be run from the hydrodynamics folder. """
    if options.daemon == True:
        evolve_plummer(daemon=True)
        smash_plummers(daemon=True)
        create_animation(daemon=True)

    else:
        evolve_plummer()
        smash_plummers()
        create_animation()
    
def evolve_plummer(daemon=False):
    if daemon == True:
#        command = 'amuse evolve.py -N 2000 -n 300 -t 0.5 -H evolveN2000n300t0.5.hdf5 -B bodyN2000n300.hdf5'
#        subprocess.call(command, shell=True)
        command = 'amuse evolve.py -N 1000 -n 100 -t 0.5 -H evolveN1000n100t0.5.hdf5 -B bodyN1000n100.hdf5'
        subprocess.call(command, shell=True)
    else:
        command = 'mpiexec amuse evolve.py -N 2000 -n 300 -t 0.5 -H evolveN2000n300t0.5.hdf5 -B bodyN2000n300.hdf5'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse evolve.py -N 1000 -n 100 -t 0.5 -H evolveN1000n100t1.hdf5 -B bodyN1000n100.hdf5'
        subprocess.call(command, shell=True)

def smash_plummers(daemon=False):
    if daemon == True:
#        command = 'amuse evolve.py -n 300 -t 0.5 -H no_velocity.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5'
#        subprocess.call(command, shell=True)
#        command = 'amuse evolve.py -n 200 -t 1.5 -H smash_vx20.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 20'
#        subprocess.call(command, shell=True)
#        command = 'amuse evolve.py -n 200 -t 0.5 1 -H smash_vx76.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 76'
#        subprocess.call(command, shell=True)
#        command = 'amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 5'
#        subprocess.call(command, shell=True)
        command = 'amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy10 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 10'
        subprocess.call(command, shell=True)
    else:
        command = 'mpiexec amuse evolve.py -n 300 -t 0.5 -H no_velocity.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse evolve.py -n 200 -t 1.5 -H smash_vx20.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 20'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse evolve.py -n 200 -t 0.5 1 -H smash_vx76.hdf5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 76'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy5 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 5'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy10 -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 10'
        subprocess.call(command, shell=True)

def create_animation(daemon=False):
    if daemon == True:
#        command = 'amuse animation.py -f hydroresults/evolveN2000n300t0.5.hdf5 -r 2'
#        subprocess.call(command, shell=True)
#        command = 'amuse animation.py -f hydroresults/no_velocity.hdf5 -r 3'
#        subprocess.call(command, shell=True)
#        command = 'amuse animation.py -f hydroresults/smash_vx20.hdf5 -r 4'
#        subprocess.call(command, shell=True)
#        command = 'amuse animation.py -f hydroresults/smash_vx76.hdf5 -r 8'
#        subprocess.call(command, shell=True)
#        command = 'amuse animation.py -f hydroresults/smash_vx10vy5 -r 4'
#        subprocess.call(command, shell=True)
        command = 'amuse animation.py -f hydroresults/smash_vx10vy10 -r 4'
        subprocess.call(command, shell=True)
    else:
        command = 'mpiexec amuse animation.py -f hydroresults/evolveN2000n300t0.5.hdf5 -r 2'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse animation.py -f hydroresults/no_velocity.hdf5 -r 3'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse animation.py -f hydroresults/smash_vx20.hdf5 -r 4'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse animation.py -f hydroresults/smash_vx76.hdf5 -r 8'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse animation.py -f hydroresults/smash_vx10vy5 -r 4'
        subprocess.call(command, shell=True)
        command = 'mpiexec amuse animation.py -f hydroresults/smash_vx10vy10 -r 4'
        subprocess.call(command, shell=True)

def parse_sysargs(sysargs):
    parser = OptionParser()

    parser.add_option("-d", "--daemon", action="store_true",\
                      dest="daemon", default=False,\
                      help="Flag this if MPI daemon is running.")

    options, args = parser.parse_args(sysargs[1:])
    return options

if __name__ == "__main__":
    options = parse_sysargs(sys.argv)
    main(options)

