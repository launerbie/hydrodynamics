

####################   Evolve plummer models like this: #######################
# % mpiexec amuse hydro_fi.py -N 4000 -n 100 -t 1.0 -H evolveN4000n100t1.hdf5
# % mpiexec amuse hydro_fi.py -N 2000 -n 300 -t 0.5 -H evolveN2000n300t0.5.hdf5
# % mpiexec amuse hydro_fi.py -N 1000 -n 100 -t 0.5 -H evolveN1000n100t1.hdf5
#
# This will generate 3 bodies in the directory 'bodies/':
# bodies/bodyN4000n100.hdf5
# bodies/bodyN2000n300.hdf5
# bodies/bodyN1000n100.hdf5
#
# And 3 hydroresults in the directory 'hydroresults/'
# hydroresults/evolveN4000n100t1.hdf5 
# hydroresults/evolveN2000n300t0.5.hdf5 
# hydroresults/evolveN1000n100t1.hdf5
# 
##############################################################################

###    Smash plummer models like this: #####
#
# % mpiexec amuse hydro_fi.py -n 200 -t 0.5 -H smash_vx20.hdf5 \
#   -p body_N1000n50_1.hdf5 -q body_N1000n50_2.hdf5 --vx 20 --vy 0 --vz 0 
#
# This will 
#
#mpiexec amuse64 animation.py -f hydroresults/evolveN2000n300t0.5.hdf5 -r 8
#
#


#mpiexec amuse hydro_fi.py -N 1000 -n 100 -t 0.5 -B testbody1.hdf5

mpiexec amuse64 hydro_fi.py -n 100 -t 1.5 -H smash_vx15vy5.hdf5 \
-p bodies/testbody1.hdf5 -q bodies/testbody1.hdf5 --vx 15 --vy 5 --vz 0 

mpiexec amuse64 animation.py -f hydroresults/smash_vx15vy5.hdf5 -r 4



