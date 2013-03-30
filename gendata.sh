
####################   Evolve plummer models like this: #######################
# % mpiexec amuse evolve.py -N 4000 -n 100 -t 1.0 -H evolveN4000n100t1.hdf5 
# % mpiexec amuse evolve.py -N 2000 -n 300 -t 0.5 -H evolveN2000n300t0.5.hdf5 -B bodyN2000n300.hdf5
# % mpiexec amuse evolve.py -N 1000 -n 100 -t 0.5 -B bodyN1000n100.hdf5
#
# This will generate the following bodies in the  default directory 'bodies/':
# bodies/bodyN2000n300.hdf5
# bodies/bodyN1000n100.hdf5
#
# And the following hydroresults in the default directory 'hydroresults/':
# hydroresults/evolveN4000n100t1.hdf5 
# hydroresults/evolveN2000n300t0.5.hdf5 
# 
##############################################################################

#mpiexec amuse evolve.py -N 2000 -n 300 -t 0.5 -H evolveN2000n300t0.5.hdf5 -B bodyN2000n300.hdf5
#mpiexec amuse evolve.py -N 1000 -n 100 -t 0.5 -H evolveN1000n100t1.hdf5 -B bodyN1000n100.hdf5



####################   Smash plummer models like this: #########################
#
# % mpiexec amuse evolve.py -n 200 -t 1 -H smash_vx20.hdf5 \
#   -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 20 --vy 0 --vz 0 
#
# % mpiexec amuse evolve.py -n 200 -t 1 -H smash_vx76.hdf5 \
#   -p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 76 --vy 0 --vz 0 
#
# This will generate the following hydroresults in the default directory 'hydroresults/'
# hydroresults/smash_vx20.hdf5 
# hydroresults/smash_vx76.hdf5
#
#
##################################################################################

mpiexec amuse evolve.py -n 300 -t 0.5 -H no_velocity.hdf5 \
-p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 

#mpiexec amuse evolve.py -n 200 -t 1.5 -H smash_vx20.hdf5 \
#-p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 20  
#
#mpiexec amuse evolve.py -n 200 -t 0.5 1 -H smash_vx76.hdf5 \
#-p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 76  
#
#mpiexec amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy5 \
#-p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 5  
#
#mpiexec amuse evolve.py -n 200 -t 1.5 1 -H smash_vx10vy10 \
#-p bodies/bodyN1000n100.hdf5 -q bodies/bodyN1000n100.hdf5 --vx 10 --vy 10  



####################   Create animations  ########################################
#
# % mpiexec amuse animation.py -f hydroresults/evolveN2000n300t0.5.hdf5 -r 2
# % mpiexec amuse animation.py -f hydroresults/smash_vx20.hdf5 -r 4
# % mpiexec amuse animation.py -f hydroresults/smash_vx76.hdf5 -r 8
#
# This will generate the following animations in the directory 'movies':
#
# evolveN2000n300t0.5r2.0.mp4
# smash_vx20r4.0.mp4 
# smash_vx76r8.0.mp4 
#
#################################################################################

#mpiexec amuse animation.py -f hydroresults/evolveN2000n300t0.5.hdf5 -r 2
mpiexec amuse animation.py -f hydroresults/no_velocity.hdf5 -r 3
#mpiexec amuse animation.py -f hydroresults/smash_vx20.hdf5 -r 4
#mpiexec amuse animation.py -f hydroresults/smash_vx76.hdf5 -r 8
#mpiexec amuse animation.py -f hydroresults/smash_vx10vy5 -r 4
#mpiexec amuse animation.py -f hydroresults/smash_vx10vy10 -r 4




